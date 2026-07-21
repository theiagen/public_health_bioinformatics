/* ============================================================================
 * PRIMER for non-JavaScript developers (shared by all table-*.js files)
 * ============================================================================
 *
 * These scripts run in the reader's browser and enhance the plain HTML tables
 * that Markdown/Zensical produce — adding search, sorting, column resizing, and
 * a scroll-to-top button.
 *
 * HOW THE FILES COOPERATE
 *   This file (table-utils.js) loads FIRST and defines a single shared object,
 *   window.mdTables, holding the helpers and the one constant that more than one
 *   feature needs. Every other table script reads from mdTables, so there is one
 *   copy of each shared helper and one page-load registration style:
 *     - table-sort.js       click-to-sort (via the Tablesort library)
 *     - table-resize.js     drag-to-resize columns
 *     - table-scrolltop.js  floating "scroll to top of table" button
 *     - table-search.js     search box + reset button
 *   Each feature also publishes its own reset + "is this active?" check onto
 *   mdTables (e.g. resetSort/isSortActive), so table-search.js's reset button can
 *   clear all of them without reaching into their internals.
 *
 *   The scripts avoid stepping on each other:
 *     - none of them add a `class` to the <table> (they use data- attributes),
 *       so the `table:not([class])` selector and the CSS still match, and
 *     - resize handle drags call stopPropagation so a resize never triggers a
 *       sort.
 * ============================================================================ */

/*
 * Shared foundation for the table scripts. Defines window.mdTables; must be
 * listed before the other table-*.js files in mkdocs.yml (extra_javascript).
 */
(function () {
  // px — the minimum width a column may be dragged to (matches the [data-manual]
  // min-width in extra.css). Shared by table-resize.js (the drag floor) and used
  // when computing a table's intrinsic minimum width.
  const MIN_COLUMN_WIDTH = 60;

  // Register fn to run on first page load AND after every Material "instant
  // navigation" page change, falling back to DOMContentLoaded when the document$
  // observable is unavailable. Every table script registers through this so the
  // lifecycle wiring lives in exactly one place.
  function onPageLoad(fn) {
    if (typeof document$ !== "undefined" && document$.subscribe) {
      document$.subscribe(fn);
    } else {
      document.addEventListener("DOMContentLoaded", fn);
    }
  }

  // The header row of a table: the first row of <thead> if present, else the
  // table's first row. Returns null for an empty table.
  function getHeaderRow(table) {
    return (table.tHead && table.tHead.rows[0]) || table.rows[0] || null;
  }

  // Every plain content table we enhance — the same set tablesort targets. The
  // `:not([class])` guard skips tables that Markdown authors have given an
  // explicit class (those opt out of these enhancements).
  function getContentTables() {
    return document.querySelectorAll("article table:not([class])");
  }

  // Idempotently wrap a table in a `.md-table-wrap` scroll container and return
  // that wrapper. Both table-resize.js and table-scrolltop.js need the wrapper;
  // whichever runs first creates it and the other reuses it.
  function getTableWrap(table) {
    let wrap = table.parentElement;
    if (!wrap || !wrap.classList.contains("md-table-wrap")) {
      wrap = document.createElement("div");
      wrap.className = "md-table-wrap";
      table.parentNode.insertBefore(wrap, table);
      wrap.appendChild(table);
    }
    return wrap;
  }

  window.mdTables = {
    MIN_COLUMN_WIDTH,
    onPageLoad,
    getHeaderRow,
    getContentTables,
    getTableWrap,
  };
})();
