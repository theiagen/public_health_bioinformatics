/*
 * Makes each plain content table sortable: clicking a header sorts by that
 * column. The sorting itself is done by the Tablesort library (loaded
 * separately); this file just attaches it to the right tables on every page.
 *
 * Depends on window.mdTables (table-utils.js) for the page-load hook and the
 * content-table selector, and on the Tablesort library. It publishes
 * mdTables.resetSort(table) and mdTables.isSortActive(table) so table-search.js's
 * reset button can clear the sort without knowing how Tablesort marks it.
 *
 * Tables with a single data row are skipped — there is nothing to sort.
 */
(function () {
  const { onPageLoad, getContentTables } = window.mdTables;

  // Tablesort marks the sorted column by setting aria-sort on its header; both
  // helpers below key off that single marker.
  function isSortActive(table) {
    return !!table.querySelector("th[aria-sort]");
  }

  // Clear Tablesort's markers so the next header click sorts fresh (ascending).
  // Note: this only removes the sort indicator; it does not restore the original
  // row order — table-search.js's reset button handles that.
  function resetSort(table) {
    table.querySelectorAll("th[aria-sort]").forEach((th) => {
      th.removeAttribute("aria-sort");
    });
  }

  function enableSorting() {
    getContentTables().forEach(function (table) {
      // Skip tables with a single data row — there is nothing to sort, so the
      // sort affordance (the visual cue that a header is clickable to sort) would
      // be pointless (e.g. Quick Facts tables).
      const body = table.tBodies[0];
      if (!body || body.rows.length <= 1) return;
      new Tablesort(table);
    });
  }

  window.mdTables.resetSort = resetSort;
  window.mdTables.isSortActive = isSortActive;

  onPageLoad(enableSorting);
})();
