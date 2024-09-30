function loadDynamicTable(csvPath, tableElementId, enablePostbackSafe = true) {
    Papa.parse(csvPath, {
      download: true,
      header: true,
      complete: function(results) {
        var data = results.data;
        var headers = Object.keys(data[0]);  
  
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
          $(`#${tableElementId}`).DataTable().clear().destroy();
        }
  
        // Initialize DataTable
        var table = $(`#${tableElementId}`).DataTable({
          data: tableRows,
          columns: headers.map(function(header) {
            return { title: header };
          }),
          paging: false,
          searching: true,
          ordering: true,
          order: [],
          autoWidth: false,
          scrollY: "600px",  // Scrollable height
          scrollX: true,      // Enable horizontal scroll if needed
          scrollCollapse: true,  // Collapse the scroll area when fewer rows
          responsive: true,
          stateSave: true
        });
  
        // Adjust columns after initialization
        table.columns.adjust();
  
        // Initialize colResizable after DataTable
        $(`#${tableElementId}`).colResizable({
          liveDrag: true,
          headerOnly: false,
          postbackSafe: enablePostbackSafe
        });
      }
    });
  }
  