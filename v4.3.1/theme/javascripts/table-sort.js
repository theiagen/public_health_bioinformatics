/*
 * Make tables sortable with Tablesort
 */
(function () {
  // functions from table-utils.js
  const { onPageLoad, getContentTables } = window.mdTables;

  // Tablesort marks a sorted column by setting aria-sort on its header
  function isSortActive(table) {
    return !!table.querySelector("th[aria-sort]");
  }

  // clear any Tablesort's markers so the sort tag is on the sorted column
  function resetSort(table) {
    table.querySelectorAll("th[aria-sort]").forEach((th) => {
      th.removeAttribute("aria-sort");
    });
  }

  // set up sorting on all applicable tables
  function enableSorting() {
    getContentTables().forEach(function (table) {
      // add Tablesort to all tables (except ones with a single row)
      const body = table.tBodies[0];
      if (!body || body.rows.length <= 1) return;
      new Tablesort(table);
    });
  }

  // enable the reset button to reset sorting
  window.mdTables.resetSort = resetSort;
  window.mdTables.isSortActive = isSortActive;

  // run this when the page loads
  onPageLoad(enableSorting);
})();
