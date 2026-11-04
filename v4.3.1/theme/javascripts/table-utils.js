/*
 * shared utilities for the other javascripts;
 * must be listed before the other files in mkdocs.yml (extra_javascript).
 */
(function () {
  // the minimum width a column can be dragged to
  const MIN_COLUMN_WIDTH = 60;

  // run a function whenever the page is loaded
  function onPageLoad(fn) {
    if (typeof document$ !== "undefined" && document$.subscribe) {
      document$.subscribe(fn);
    } else {
      document.addEventListener("DOMContentLoaded", fn);
    }
  }

  // get the header row of a table: the first row of <thead> if present, else the
  // table's first row. return null for an empty table
  function getHeaderRow(table) {
    return (table.tHead && table.tHead.rows[0]) || table.rows[0] || null;
  }

  // grab every plain content table
  // `:not([class])` skips tables that we've given an explicit class
  function getContentTables() {
    return document.querySelectorAll("article table:not([class])");
  }

  // wrap a table in a `.md-table-wrap` scroll container and return it
  function getTableWrap(table) {
    let wrap = table.parentElement;
    // only add the scroll container if it doesn't have it yet
    if (!wrap || !wrap.classList.contains("md-table-wrap")) {
      wrap = document.createElement("div");
      wrap.className = "md-table-wrap";
      table.parentNode.insertBefore(wrap, table);
      wrap.appendChild(table);
    }
    return wrap;
  }

  // make these functions accessible by other javascripts
  window.mdTables = {
    MIN_COLUMN_WIDTH,
    onPageLoad,
    getHeaderRow,
    getContentTables,
    getTableWrap,
  };
})();
