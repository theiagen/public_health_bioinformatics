document$.subscribe(function() {
  var tables = document.querySelectorAll("article table:not([class])")
  tables.forEach(function(table) {
    // Skip tables with a single data row — there is nothing to sort, so the
    // sort affordance would be pointless (e.g. Quick Facts tables).
    var body = table.tBodies[0]
    if (!body || body.rows.length <= 1) return
    new Tablesort(table)
  })
})