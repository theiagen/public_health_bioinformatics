import yaml
from cryptography.fernet import Fernet
from js import Blob, Object, URL, navigator
from pyodide.ffi import to_js
from pyscript import document, when


def value(name):
    """Strip all input fields of leading and trailing whitespace and report None if blank."""
    return document.getElementById(f"cb-{name}").value.strip() or None


def checked(name):
    """Reads the checked state of a checkbox on the page and returns it as a boolean value"""
    return document.getElementById(f"cb-{name}").checked


def build_config(encrypted_password):
    """
    Mimicking the SeqSender submission wizard config file builder in: "https://cdcgov.github.io/seqsender/".
    Some fields left out for better compatibility with Terra workflow implementation.
    """
    biosample = checked("biosample")
    genbank = checked("genbank")
    ncbi = {
        "Username": value("username"),
        "Password": encrypted_password,
        "Spuid_Namespace": value("spuid_namespace"),
    }
    if biosample:
        ncbi["BioSample_Package"] = value("biosample_package")

    # Intentionally modified/omitted the following fields (seqsender treats a missing value as off):
    # - GenBank_Auto_Remove_Failed_Samples
    # - Specified_Release_Date
    # - Link_Sample_Between_NCBI_Databases
    # - Add_Definition_Line_Accessions
    ncbi.update({
        "Publication_Title": value("publication_title") if genbank else None,
        "Publication_Status": value("publication_status") if genbank else None,
        "Specified_Release_Date": None,  # requires the key; empty releases the data immediately once NCBI processes it
        "Link_Sample_Between_NCBI_Databases": checked("link_samples") if genbank else True, # only GenBank uses this; without GenBank it's always True, since seqsender always links BioSample and SRA
        "Add_Definition_Line_Accessions": checked("definition_line_accessions") if genbank else None,
        "Description": {
            "Organization": {
                "Role": value("role"),
                "Type": value("type"),
                "Name": value("organization"),
                "Address": {
                    "Affil": value("affil"),
                    "Div": value("div"),
                    "Street": value("street"),
                    "City": value("city"),
                    "Sub": value("sub"),
                    "Postal_Code": value("postal_code"),
                    "Country": value("country"),
                    "Email": value("email"),
                    "Phone": value("phone"),
                },
                "Submitter": {
                    "Email": value("submitter_email"),
                    "Alt_Email": value("submitter_alt_email"),
                    "Name": {
                        "First": value("first_name"),
                        "Last": value("last_name")
                    },
                },
            },
        },
    })
    return {"Submission": {"NCBI": ncbi}}


@when("change", ".cb-database")
def toggle_database_fields(event):
    # a disabled fieldset is hidden (extra.css) and skipped by the browser's required-field checks
    for database in ("biosample", "genbank"):
        document.getElementById(f"cb-{database}-fields").disabled = not checked(database)


@when("submit", "#cb-form")
def create_config(event):
    event.preventDefault()
    document.getElementById("cb-error").hidden = True
    document.getElementById("cb-result").hidden = True

    # passwords aren't trimmed; spaces can be part of them
    password = document.getElementById("cb-password").value

    # reusing a key lets one key unlock several config files
    key = value("key") or Fernet.generate_key().decode()
    try:
        # key gets baked into the Fernet object
        fernet = Fernet(key)
    except ValueError:
        return show_error('Invalid SeqSender key. A key is 44 characters long and ends with "=".')

    # like SeqSender's `load_credentials`: the encrypted password is bytes, which yaml.safe_dump writes as !!binary
    config_yaml = yaml.safe_dump(build_config(fernet.encrypt(password.encode())), sort_keys=False)
    blob = Blob.new(to_js([config_yaml]), to_js({"type": "text/yaml"}, dict_converter=Object.fromEntries))
    document.getElementById("cb-download").href = URL.createObjectURL(blob)
    document.getElementById("cb-key-output").textContent = key
    document.getElementById("cb-copy-key").textContent = "Copy key"
    document.getElementById("cb-result").hidden = False


@when("click", "#cb-copy-key")
def copy_key(event):
    navigator.clipboard.writeText(document.getElementById("cb-key-output").textContent)
    event.target.textContent = "Copied"


def show_error(message):
    error = document.getElementById("cb-error")
    error.lastElementChild.textContent = message
    error.hidden = False


# the Create button stays disabled until Python has loaded
document.getElementById("cb-loading").hidden = True
document.getElementById("cb-create").disabled = False
