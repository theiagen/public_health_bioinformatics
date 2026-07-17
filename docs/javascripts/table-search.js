/*
 * Adds a search box + reset button above each ".searchable-table" and filters
 * its rows as the reader types.
 *
 * Terms used below:
 *   - DOM (Document Object Model): the browser's live, in-memory tree of the
 *     page. JavaScript reads and changes the page by walking this tree.
 *   - document$ : a Material for MkDocs "observable" (a signal you can subscribe
 *     to). document$.subscribe(fn) runs fn once when the page loads AND again
 *     after every "instant navigation" — Material swaps page content without a
 *     full reload, so an ordinary page-load event would not fire on later pages.
 *   - MutationObserver: a browser tool that watches an element for changes (here,
 *     an attribute being set) and calls a function whenever they happen.
 *   - requestAnimationFrame(fn): asks the browser to run fn just before the next
 *     screen repaint (when the browser redraws the page). Used to batch rapid
 *     work — e.g. filtering while the reader types — into one update per frame
 *     instead of one per keystroke.
 *   - Function "hoisting": a `function foo() {}` can be called from code written
 *     above it, so helpers are defined in a readable order, not strictly before
 *     first use.
 *
 * What this file does specifically:
 *   - caches each row's text once, so typing never re-reads the DOM;
 *   - runs the filter on the next animation frame, so fast keystrokes are
 *     batched into a single update;
 *   - narrows incrementally: when the query grows (a character is typed), only
 *     rows that are still visible are re-checked, because a row already hidden
 *     cannot match a longer query;
 *   - watches the sorted-column marker (the aria-sort attribute) with a
 *     MutationObserver so the reset button appears the moment a column is
 *     sorted;
 *   - re-applies itself on every page via document$.
 *
 * (A fuller primer on these shared concepts also lives at the top of
 * table-resize.js.)
 */

// Reset icon (circular arrow) used by the per-table reset button.
const RESET_ICON =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
    '<path d="M12 5V1L7 6l5 5V7c3.31 0 6 2.69 6 6s-2.69 6-6 6-6-2.69-6-6H4c0 4.42 3.58 8 8 8s8-3.58 8-8-3.58-8-8-8z"/>' +
    '</svg>';

function addTableSearch() {
    const containers = document.querySelectorAll('.searchable-table');

    containers.forEach((container) => {
        const table = container.querySelector('table');
        if (!table) {
            console.log('Table not found within container.');
            return;
        }

        // Don't add a second toolbar if this container is already enhanced.
        if (container.querySelector('.table-toolbar')) {
            return;
        }

        // ----- Toolbar: search box (left) + reset button (right) -----
        const toolbar = document.createElement('div');
        toolbar.classList.add('table-toolbar');

        const searchInput = document.createElement('input');
        searchInput.setAttribute('type', 'search');
        searchInput.setAttribute('placeholder', 'Search');
        searchInput.classList.add('table-search-input');

        // Reset button: clears the search, sort order, and manual column widths.
        // Pinned to the far right via CSS and hidden until something is active.
        const resetButton = document.createElement('button');
        resetButton.type = 'button';
        resetButton.classList.add('table-reset-button', 'is-hidden');
        resetButton.title = 'Reset table';
        resetButton.setAttribute('aria-label', 'Reset table');
        resetButton.innerHTML = RESET_ICON;

        toolbar.appendChild(searchInput);
        toolbar.appendChild(resetButton);
        container.insertBefore(toolbar, container.firstChild);

        // ----- Row bookkeeping -----
        // Capture the original row order so reset can restore it (tablesort
        // reorders these same <tr> nodes in place, so the references stay valid).
        const tbody = table.tBodies[0];
        const originalOrder = tbody ? Array.from(tbody.rows) : [];

        // Live collection of the data rows (header lives in <thead>, so this is
        // body-only). Reused across keystrokes instead of re-querying each time.
        const bodyRows = tbody ? tbody.rows : [];

        // Return a row's uppercased searchable text, caching it on the element.
        // textContent (not innerText) avoids forcing a layout on every read, and
        // the cache means repeated keystrokes never touch the DOM for reading.
        // The cache rides on the <tr>, so tablesort/reset reordering keeps it.
        function rowSearchText(row) {
            if (row._searchText === undefined) {
                row._searchText = row.textContent.toUpperCase();
            }
            return row._searchText;
        }

        function showAllRows() {
            for (let i = 0; i < bodyRows.length; i++) {
                bodyRows[i].style.display = '';
            }
        }

        // ----- Reset-button visibility -----
        function isSortActive() {
            // tablesort marks the sorted column with the aria-sort attribute
            return !!table.querySelector('th[aria-sort]');
        }

        function isResizeActive() {
            return table.dataset.manual === 'true';
        }

        // Show the reset button only once a search, sort, or resize is active.
        function updateResetVisibility() {
            const active =
                searchInput.value.trim() !== '' ||
                isSortActive() ||
                isResizeActive();
            resetButton.classList.toggle('is-hidden', !active);
        }

        // ----- Search filtering -----
        // State for incremental narrowing (see applyFilter).
        let filterScheduled = false;
        let lastFilter = '';
        let visibleRows = null; // null = unknown, must re-scan every row

        // Apply the current filter against the cached per-row text. Only writes
        // style.display when a row's visibility actually flips, so unchanged rows
        // never invalidate layout.
        //
        // Incremental narrowing: when the new query CONTAINS the previous one
        // (e.g. the user typed another character), any row already hidden cannot
        // match the longer query, so we only re-check the rows still visible.
        // When the query widens (characters removed) we re-scan every row.
        function applyFilter() {
            filterScheduled = false;
            const filter = searchInput.value.toUpperCase();

            const narrowing = visibleRows !== null && filter.includes(lastFilter);
            const rows = narrowing ? visibleRows : bodyRows;
            const nextVisible = [];

            for (let i = 0; i < rows.length; i++) {
                const row = rows[i];
                const match = filter === '' || rowSearchText(row).includes(filter);

                if (match) {
                    if (row.style.display === 'none') row.style.display = '';
                    nextVisible.push(row);
                } else if (row.style.display !== 'none') {
                    row.style.display = 'none';
                }
            }

            visibleRows = nextVisible;
            lastFilter = filter;
            updateResetVisibility();
        }

        // Filter as the user types. Rapid keystrokes coalesce into a single pass
        // on the next animation frame instead of one full pass per keypress.
        searchInput.addEventListener('input', function () {
            if (filterScheduled) return;
            filterScheduled = true;
            requestAnimationFrame(applyFilter);
        });

        // ----- React to sort / resize so the reset button appears -----
        // Tablesort sets the aria-sort attribute on the sorted header; watch for
        // that so the reset button appears as soon as a column is sorted.
        const thead = table.tHead || table;
        const sortObserver = new MutationObserver(updateResetVisibility);
        sortObserver.observe(thead, {
            subtree: true,
            attributes: true,
            attributeFilter: ['aria-sort'],
        });
        document.addEventListener('table-resized', updateResetVisibility);

        // ----- Reset: clear search, restore order, drop sort + resize -----
        resetButton.addEventListener('click', function () {
            searchInput.value = '';

            if (tbody) {
                originalOrder.forEach((row) => tbody.appendChild(row));
            }

            // Clear tablesort's aria-sort markers so the next header click sorts
            // fresh (ascending).
            table.querySelectorAll('th[aria-sort]').forEach((th) => {
                th.removeAttribute('aria-sort');
            });

            showAllRows();
            // Reset narrowing state so the next filter re-scans every row.
            visibleRows = null;
            lastFilter = '';
            resetTableResize(table);
            updateResetVisibility();
        });
    });
}

// Run on load and on every instant-navigation page change; addTableSearch is
// idempotent per container. Uses the same hook as table-resize.js/table-sort.js.
if (typeof document$ !== "undefined" && document$.subscribe) {
    document$.subscribe(addTableSearch);
} else {
    document.addEventListener("DOMContentLoaded", addTableSearch);
}
