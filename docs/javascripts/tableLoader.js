function loadSimpleTable(csvPath, tableElementId) {
  console.log("Attempting to load table with CSV:", csvPath);
  
  Papa.parse(csvPath, {
    download: true,
    header: true,
    complete: function(results) {
      console.log("CSV data fetched successfully:", results);
      
      var data = results.data;
      if (data.length === 0) {
        console.error("No data found in the CSV");
        return;
      }

      var headers = Object.keys(data[0]);
      console.log("CSV headers:", headers);

      // Build table headers
      var headerHTML = '<tr>';
      headers.forEach(function(header) {
        headerHTML += '<th>' + header + '</th>';
      });
      headerHTML += '</tr>';
      document.querySelector(`#${tableElementId} thead`).innerHTML = headerHTML;

      // Build table data from CSV
      var tableRows = [];
      data.forEach(function(row) {
        var rowArray = [];
        headers.forEach(function(header) {
          var cellContent = row[header] ? row[header] : '';
          rowArray.push(cellContent);
        });
        tableRows.push(rowArray);
      });

      // Ensure DataTable is destroyed and reinitialized if it exists
      if ($.fn.DataTable.isDataTable(`#${tableElementId}`)) {
        console.log("Reinitializing DataTable");
        $(`#${tableElementId}`).DataTable().clear().destroy();
      }

      // Initialize basic DataTable
      $(`#${tableElementId}`).DataTable({
        data: tableRows,
        columns: headers.map(function(header) {
          return { title: header };
        }),
        paging: false,
        searching: true,
        ordering: true,
        autoWidth: false,  // Disable autoWidth to control column sizes
        scrollX: true,     // Enable horizontal scrolling
        stateSave: true    // Save state like column widths
      });
    }
  });
}

document.addEventListener('DOMContentLoaded', function() {
  // Initialize the simple table when the DOM is ready
  loadSimpleTable('/public_health_bioinformatics/data/assembly_fetch_input.csv', 'inputTable');
});
