from datetime import datetime
def on_config(config, **kwargs):
    config.copyright = f"© 2022-{datetime.now().year} <a href=\"https://www.theiagen.com\" target=\"_blank\" rel=\"noopener\">Theiagen Genomics</a>"