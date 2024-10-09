function addTableSearch() {
    // Array of table IDs to add search functionality to
    const tableIds = ['inputs-table', 'outputs-table'];

    tableIds.forEach((tableId) => {
        // Find the element with the given ID
        const container = document.getElementById(tableId);

        if (container) {
            // Find the table within this container
            const table = container.querySelector('table');

            if (table) {
                // Ensure we don't add multiple search boxes
                if (!container.querySelector('input[type="search"]')) {
                    // Create the search input element
                    const searchInput = document.createElement("input");
                    searchInput.setAttribute("type", "search");
                    searchInput.setAttribute("placeholder", "Search table...");
                    searchInput.classList.add('table-search-input');
                    searchInput.style.marginBottom = "10px";
                    searchInput.style.display = "block";

                    // Insert the search input before the table
                    container.insertBefore(searchInput, container.firstChild);

                    // Add event listener for table search
                    searchInput.addEventListener("input", function () {
                        const filter = searchInput.value.toUpperCase();
                        const rows = table.getElementsByTagName("tr");

                        for (let i = 1; i < rows.length; i++) { // Skip header row
                            const cells = rows[i].getElementsByTagName("td");
                            let match = false;

                            for (let j = 0; j < cells.length; j++) {
                                if (cells[j].innerText.toUpperCase().includes(filter)) {
                                    match = true;
                                    break;
                                }
                            }

                            rows[i].style.display = match ? "" : "none";
                        }
                    });
                }
            } else {
                console.log('Table not found within element with ID "' + tableId + '".');
            }
        } else {
            console.log('Element with ID "' + tableId + '" not found.');
        }
    });
}

// Run on page load
addTableSearch();

// Reapply search bar on page change
function observeDOMChanges() {
    const targetNode = document.querySelector('body');
    const config = { childList: true, subtree: true };

    const observer = new MutationObserver(() => {
        addTableSearch();
    });

    observer.observe(targetNode, config);
}

observeDOMChanges();
