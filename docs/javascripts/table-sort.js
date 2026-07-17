/*
 * Makes each plain content table sortable: clicking a header sorts by that
 * column. The sorting itself is done by the Tablesort library (loaded
 * separately); this file just attaches it to the right tables on every page.
 *
 * Terms used below:
 *   - DOM (Document Object Model): the browser's live, in-memory tree of the
 *     page. JavaScript reads and changes the page by walking this tree.
 *   - document.querySelectorAll(selector): finds ALL elements in the DOM that
 *     match a CSS selector, returning a list you can loop over with forEach.
 *     (querySelector — no "All" — returns just the first match.)
 *   - document$ : a Material for MkDocs "observable" (a signal you can subscribe
 *     to). document$.subscribe(fn) runs fn once when the page loads AND again
 *     after every "instant navigation" — Material swaps page content without a
 *     full reload, so an ordinary page-load event would not fire on later pages.
 *
 * Tables with a single data row are skipped — there is nothing to sort. (A
 * fuller primer on these concepts lives at the top of table-resize.js.)
 */
document$.subscribe(function() {
  var tables = document.querySelectorAll("article table:not([class])")
  tables.forEach(function(table) {
    // Skip tables with a single data row — there is nothing to sort, so the
    // sort affordance (the visual cue that a header is clickable to sort) would
    // be pointless (e.g. Quick Facts tables).
    var body = table.tBodies[0]
    if (!body || body.rows.length <= 1) return
    new Tablesort(table)
  })
})