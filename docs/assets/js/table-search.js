function addTableSearch() {
  const tables = document.querySelectorAll("table");

  tables.forEach((table) => {
      if (!table.parentElement.querySelector('input[type="search"]')) {
          // Create the search input element
          const searchInput = document.createElement("input");
          searchInput.setAttribute("type", "search");
          searchInput.setAttribute("placeholder", "Search table..."); // Updated placeholder text
          searchInput.style.marginBottom = "10px";
          searchInput.style.display = "block"; // Ensure the input box appears above the table

          // Insert the search input before the table
          table.parentElement.insertBefore(searchInput, table);

          // Add event listener for table search
          searchInput.addEventListener("input", function () {
              const filter = searchInput.value.toUpperCase();
              const rows = table.getElementsByTagName("tr");

              // Loop through all table rows, except the first (header row)
              for (let i = 1; i < rows.length; i++) {
                  const cells = rows[i].getElementsByTagName("td");
                  let match = false;

                  // Loop through all table cells to see if there's a match
                  for (let j = 0; j < cells.length; j++) {
                      if (cells[j].innerText.toUpperCase().includes(filter)) {
                          match = true;
                          break;
                      }
                  }

                  // Show or hide rows based on matching content
                  rows[i].style.display = match ? "" : "none";
              }
          });
      }
  });
}

// Function to monitor DOM changes and reinitialize the search box when a new page is loaded
function observeDOMChanges() {
  const targetNode = document.querySelector('body'); // Monitor the body for changes
  const config = { childList: true, subtree: true }; // Look for added/removed nodes

  const observer = new MutationObserver(() => {
      addTableSearch(); // Reapply search box when a new page is loaded
  });

  // Start observing the target node for configured mutations
  observer.observe(targetNode, config);
}

// Run on initial page load
addTableSearch();

// Start observing DOM changes to handle page navigation in MkDocs
observeDOMChanges();
