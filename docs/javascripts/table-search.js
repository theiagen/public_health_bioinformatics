// Reset icon (circular arrow) used by the per-table reset button.
const RESET_ICON =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
    '<path d="M12 5V1L7 6l5 5V7c3.31 0 6 2.69 6 6s-2.69 6-6 6-6-2.69-6-6H4c0 4.42 3.58 8 8 8s8-3.58 8-8-3.58-8-8-8z"/>' +
    '</svg>';

function addTableSearch() {
    // Select all containers with the class 'searchable-table'
    const containers = document.querySelectorAll('.searchable-table');

    containers.forEach((container) => {
        // Find the table within this container
        const table = container.querySelector('table');

        if (!table) {
            console.log('Table not found within container.');
            return;
        }

        // Ensure we don't add multiple toolbars
        if (container.querySelector('.table-toolbar')) {
            return;
        }

        // Build the toolbar: reset button + search box, on one line above the table
        const toolbar = document.createElement('div');
        toolbar.classList.add('table-toolbar');

        // Search input
        const searchInput = document.createElement("input");
        searchInput.setAttribute("type", "search");
        searchInput.setAttribute("placeholder", "Filter this table…");
        searchInput.classList.add('table-search-input');

        // Reset button: clears the search and the sort order. Anchored to the
        // far right of the toolbar and hidden until a search or sort is active.
        const resetButton = document.createElement('button');
        resetButton.type = 'button';
        resetButton.classList.add('table-reset-button', 'is-hidden');
        resetButton.title = 'Reset table';
        resetButton.setAttribute('aria-label', 'Reset table');
        resetButton.innerHTML = RESET_ICON;

        // search box on the left, reset button pinned to the right (via CSS)
        toolbar.appendChild(searchInput);
        toolbar.appendChild(resetButton);

        // Insert the toolbar before the table (or its scroll wrapper)
        container.insertBefore(toolbar, container.firstChild);

        // Capture the original row order so the reset button can restore it
        // (tablesort reorders these same <tr> nodes in place, so the references
        // stay valid).
        const tbody = table.tBodies[0];
        const originalOrder = tbody ? Array.from(tbody.rows) : [];

        function showAllRows() {
            const rows = table.getElementsByTagName("tr");
            for (let i = 1; i < rows.length; i++) { // Skip header row
                rows[i].style.display = "";
            }
        }

        function isSortActive() {
            // tablesort marks the sorted column with the aria-sort attribute
            return !!table.querySelector("th[aria-sort]");
        }

        function isResizeActive() {
            return table.dataset.manual === "true";
        }

        // Show the reset button only once a search or sort has been applied
        function updateResetVisibility() {
          const active = searchInput.value.trim() !== "" ||
            isSortActive() ||
            isResizeActive();

            resetButton.classList.toggle("is-hidden", !active);
        }

        // Filter rows as the user types
        searchInput.addEventListener("input", function () {
            const filter = searchInput.value.toUpperCase();
            const rows = table.getElementsByTagName("tr");

            for (let i = 1; i < rows.length; i++) { // Skip header row
                const cells = rows[i].getElementsByTagName("td");
                let match = false;

                for (let j = 0; j < cells.length; j++) {
                    if (cells[j].innerText.toUpperCase().includes(filter)) {
                        match = true;
                        break;
                    }
                }

                rows[i].style.display = match ? "" : "none";
            }

            updateResetVisibility();
        });

        // Tablesort sets the aria-sort attribute on the sorted header; watch for
        // that so the reset button appears as soon as a column is sorted.
        const thead = table.tHead || table;
        const sortObserver = new MutationObserver(updateResetVisibility);
        sortObserver.observe(thead, {
            subtree: true,
            attributes: true,
            attributeFilter: ["aria-sort"],
        });
        document.addEventListener("table-resized", updateResetVisibility);

        // Reset: clear the search text, restore original order, drop sort markers
        resetButton.addEventListener("click", function () {
            searchInput.value = "";

            if (tbody) {
                originalOrder.forEach((row) => tbody.appendChild(row));
            }

            // Clear tablesort's aria-sort markers so the next header click sorts
            // fresh (ascending)
            table.querySelectorAll("th[aria-sort]").forEach((th) => {
                th.removeAttribute("aria-sort");
            });

          showAllRows();
          resetTableResize(table);
          updateResetVisibility();
        });
    });
}

// Run on page load
addTableSearch();

// Reapply search bar on page change
function observeDOMChanges() {
    const targetNode = document.querySelector('body');
    const config = { childList: true, subtree: true };

    const observer = new MutationObserver(() => {
        addTableSearch();
    });

    observer.observe(targetNode, config);
}

observeDOMChanges();
