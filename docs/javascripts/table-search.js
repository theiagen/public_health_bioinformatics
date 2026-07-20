/*
 * Adds a toolbar above each ".searchable-table" containing a search box (which
 * filters the table's rows as the reader types) and a reset button (which clears
 * the search, sort, and manual column widths at once).
 *
 * This file owns TWO things, kept in clearly separated sections below:
 *   1. SEARCH — the search input and the incremental row-filtering it drives.
 *      This is the file's core responsibility.
 *   2. RESET BUTTON — the reset button and the "is anything active?" logic that
 *      shows/hides it. Reset is cross-cutting: it clears search (owned here) plus
 *      sort and resize, which it delegates to the other features via
 *      window.mdTables (resetSort/resetResize + isSortActive/isResizeActive) so
 *      this file never reaches into their internals.
 *
 * Search performance notes:
 *   - caches each row's text once, so typing never re-reads the DOM;
 *   - runs the filter on the next animation frame, so fast keystrokes are
 *     batched into a single update;
 *   - narrows incrementally: when the query grows (a character is typed), only
 *     rows that are still visible are re-checked, because a row already hidden
 *     cannot match a longer query.
 *
 * Depends on window.mdTables (table-utils.js) for the page-load hook and the
 * per-feature reset/active-check helpers.
 */
(function () {
  const { onPageLoad, resetSort, isSortActive, resetResize, isResizeActive } =
    window.mdTables;

  // Reset icon (circular arrow) used by the per-table reset button.
  const RESET_ICON =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
    '<path d="M12 5V1L7 6l5 5V7c3.31 0 6 2.69 6 6s-2.69 6-6 6-6-2.69-6-6H4c0 4.42 3.58 8 8 8s8-3.58 8-8-3.58-8-8-8z"/>' +
    '</svg>';

  function addTableSearch() {
    const containers = document.querySelectorAll(".searchable-table");

    containers.forEach((container) => {
      const table = container.querySelector("table");
      if (!table) {
        console.log("Table not found within container.");
        return;
      }

      // Don't add a second toolbar if this container is already enhanced.
      if (container.querySelector(".table-toolbar")) {
        return;
      }

      // ----- Toolbar: search box (left) + reset button (right) -----
      const toolbar = document.createElement("div");
      toolbar.classList.add("table-toolbar");

      const searchInput = document.createElement("input");
      searchInput.setAttribute("type", "search");
      searchInput.setAttribute("placeholder", "Search");
      searchInput.classList.add("table-search-input");

      // Reset button: clears the search, sort order, and manual column widths.
      // Pinned to the far right via CSS and hidden until something is active.
      const resetButton = document.createElement("button");
      resetButton.type = "button";
      resetButton.classList.add("table-reset-button", "is-hidden");
      resetButton.title = "Reset table";
      resetButton.setAttribute("aria-label", "Reset table");
      resetButton.innerHTML = RESET_ICON;

      toolbar.appendChild(searchInput);
      toolbar.appendChild(resetButton);
      container.insertBefore(toolbar, container.firstChild);

      // ----- Row bookkeeping (shared by search + reset) -----
      // Capture the original row order so reset can restore it (tablesort
      // reorders these same <tr> nodes in place, so the references stay valid).
      const tbody = table.tBodies[0];
      const originalOrder = tbody ? Array.from(tbody.rows) : [];

      // Live collection of the data rows (header lives in <thead>, so this is
      // body-only). Reused across keystrokes instead of re-querying each time.
      const bodyRows = tbody ? tbody.rows : [];

      /* ======================================================================
       * SECTION 1 — SEARCH: filter rows as the reader types
       * ==================================================================== */

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

      // State for incremental narrowing (see applyFilter).
      let filterScheduled = false;
      let lastFilter = "";
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
          const match = filter === "" || rowSearchText(row).includes(filter);

          if (match) {
            if (row.style.display === "none") row.style.display = "";
            nextVisible.push(row);
          } else if (row.style.display !== "none") {
            row.style.display = "none";
          }
        }

        visibleRows = nextVisible;
        lastFilter = filter;
        updateResetVisibility();
      }

      // Filter as the user types. Rapid keystrokes coalesce into a single pass
      // on the next animation frame instead of one full pass per keypress.
      searchInput.addEventListener("input", function () {
        if (filterScheduled) return;
        filterScheduled = true;
        requestAnimationFrame(applyFilter);
      });

      /* ======================================================================
       * SECTION 2 — RESET BUTTON: clear search + sort + resize together
       *
       * Everything below is about the reset button, not searching. It shows the
       * button whenever ANY feature is active, and on click clears this file's
       * search state while delegating sort and resize back to their own scripts.
       * ==================================================================== */

      // Show every data row again (used by reset, after clearing the filter).
      function showAllRows() {
        for (let i = 0; i < bodyRows.length; i++) {
          bodyRows[i].style.display = "";
        }
      }

      // Show the reset button only once a search, sort, or resize is active.
      // Sort/resize state is owned by the other scripts and read via mdTables.
      function updateResetVisibility() {
        const active =
          searchInput.value.trim() !== "" ||
          isSortActive(table) ||
          isResizeActive(table);
        resetButton.classList.toggle("is-hidden", !active);
      }

      // React to sort / resize (driven by the other scripts) so the button
      // appears the moment a column is sorted or a resize is committed.
      // Tablesort sets aria-sort on the sorted header; watch for that.
      const thead = table.tHead || table;
      const sortObserver = new MutationObserver(updateResetVisibility);
      sortObserver.observe(thead, {
        subtree: true,
        attributes: true,
        attributeFilter: ["aria-sort"],
      });
      document.addEventListener("table-resized", updateResetVisibility);

      // Reset: clear this file's search + row order, then delegate sort and
      // resize back to their owning scripts.
      resetButton.addEventListener("click", function () {
        // --- search (owned here) ---
        searchInput.value = "";
        if (tbody) {
          originalOrder.forEach((row) => tbody.appendChild(row));
        }
        showAllRows();
        // Reset narrowing state so the next filter re-scans every row.
        visibleRows = null;
        lastFilter = "";

        // --- delegated to the other features ---
        resetSort(table);
        resetResize(table);

        updateResetVisibility();
      });
    });
  }

  // addTableSearch is idempotent per container (guards on an existing toolbar).
  onPageLoad(addTableSearch);
})();
