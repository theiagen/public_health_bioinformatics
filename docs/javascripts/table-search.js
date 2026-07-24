/*
 * Adds a toolbar above each ".searchable-table" containing a search box and a
 * reset button to clear any search or table visual manipulation
 */
(function () {
  // functions from table-utils.js
  const { onPageLoad } = window.mdTables;

  // icon for table reset button
  const RESET_ICON =
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" aria-hidden="true">' +
    '<path d="M12 5V1L7 6l5 5V7c3.31 0 6 2.69 6 6s-2.69 6-6 6-6-2.69-6-6H4c0 4.42 3.58 8 8 8s8-3.58 8-8-3.58-8-8-8z"/>' +
    '</svg>';

  // enable table search (and reset button) to tables
  function addTableSearch() {
    // Read these here (not at file load) so it doesn't matter that
    // table-resize.js loads after this file — by the time addTableSearch runs
    // (on page load), every table-*.js file has already published to mdTables.
    const {
      resetSort,
      isSortActive,
      resetResize,
      isResizeActive,
      scrollTableToTop,
    } = window.mdTables;

    // find all searchable table tags on the page
    const containers = document.querySelectorAll(".searchable-table");

    // for each container that has a table
    containers.forEach((container) => {
      const table = container.querySelector("table");
      if (!table) {
        return;
      }

      // only add one search bar
      if (container.querySelector(".table-toolbar")) {
        return;
      }

      // create tool bar
      const toolbar = document.createElement("div");
      toolbar.classList.add("table-toolbar");

      // add search bar (left)
      const searchInput = document.createElement("input");
      searchInput.setAttribute("type", "search");
      searchInput.setAttribute("placeholder", "Search");
      searchInput.classList.add("table-search-input");

      // add reset button (right); hidden until something is active.
      const resetButton = document.createElement("button");
      resetButton.type = "button";
      resetButton.classList.add("table-reset-button", "is-hidden");
      resetButton.title = "Reset table";
      resetButton.setAttribute("aria-label", "Reset table");
      resetButton.innerHTML = RESET_ICON;

      toolbar.appendChild(searchInput);
      toolbar.appendChild(resetButton);

      // add the toolbar before the table
      container.insertBefore(toolbar, container.firstChild);

      // save the original row order so reset can restore it (tablesort
      // reorders these same <tr> nodes in place, so the references stay valid).
      const tbody = table.tBodies[0];
      const originalOrder = tbody ? Array.from(tbody.rows) : [];

      // a live list of the rows (in whatever order the sort puts it into)
      const bodyRows = tbody ? tbody.rows : [];

/* ======================================================================
* SEARCH: filter rows as the reader types
* ==================================================================== */

      // cache a row's text in uppercase to the <tr> element so it's only
      // accessed once from the DOM and won't be lost during sorting
      function rowSearchText(row) {
        if (row._searchText === undefined) {
          row._searchText = row.textContent.toUpperCase();
        }
        return row._searchText;
      }

      // vars to allow for incremental narrowing
      let filterScheduled = false; // flag to avoid event firing when typing fast
      let lastFilter = ""; // contains what the last time the search was run
      let visibleRows = null; // the list of rows that matched last time

      // perform the search filter
      function applyFilter() {
        filterScheduled = false;
        // convert the searchbar text to uppercase
        const filter = searchInput.value.toUpperCase();

        // check to see if the search is narrowing vs widening
        // (the filter contains the same text as last time and rows are visible)
        const narrowing = visibleRows !== null && filter.includes(lastFilter);
        // if narrowing, search the visible rows; if not, search all rows
        const rows = narrowing ? visibleRows : bodyRows;
        const nextVisible = [];

        // search for the matching text in the searchable rows
        for (let i = 0; i < rows.length; i++) {
          const row = rows[i];
          // does the row have matching text?
          const match = filter === "" || rowSearchText(row).includes(filter);

          if (match) {
            // if there's a match, show the row if it was hidden
            if (row.style.display === "none") row.style.display = "";
            nextVisible.push(row);
          } else if (row.style.display !== "none") {
            // if it doesn't match, hide it if it wasn't already hidden
            row.style.display = "none";
          }
        }

        // update the list of visible rows and lastFilter
        visibleRows = nextVisible;
        lastFilter = filter;
        // indicate that a filter was set
        updateResetVisibility();
      }

      // filter as the user types. wait until they're done typing to search
      searchInput.addEventListener("input", function () {
        if (filterScheduled) return;
        filterScheduled = true;
        requestAnimationFrame(applyFilter);
      });

/* ======================================================================
* RESET BUTTON: clear search, sort, scroll, and resize together
* ==================================================================== */

      // undo any row filtering
      function showAllRows() {
        for (let i = 0; i < bodyRows.length; i++) {
          bodyRows[i].style.display = "";
        }
      }

      // only show the reset button when a search, reset, or sort was performed.
      function updateResetVisibility() {
        const active =
          searchInput.value.trim() !== "" ||
          isSortActive(table) ||
          isResizeActive(table);
        resetButton.classList.toggle("is-hidden", !active);
      }

      // react to sort / resize (driven by the other scripts) so the button
      // appears the moment a column is sorted or a resized.
      const thead = table.tHead || table;
      const sortObserver = new MutationObserver(updateResetVisibility);
      sortObserver.observe(thead, {
        subtree: true,
        attributes: true,
        attributeFilter: ["aria-sort"], // tablesort flag
      });
      document.addEventListener("table-resized", updateResetVisibility);

      // listen for when the reset button is clicked to reset
      resetButton.addEventListener("click", function () {
        // reestablish the row original order
        searchInput.value = "";
        if (tbody) {
          originalOrder.forEach((row) => tbody.appendChild(row));
        }
        showAllRows();
        // reset narrowing state so the next filter re-scans every row.
        visibleRows = null;
        lastFilter = "";

        // reset sorting, resizing, and any scrolling
        resetSort(table);
        resetResize(table);
        scrollTableToTop(table);

        // re-hide the reset button
        updateResetVisibility();
      });
    });
  }

  // when the page loads, add the table search bar
  onPageLoad(addTableSearch);
})();
