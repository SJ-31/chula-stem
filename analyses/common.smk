from datetime import date

TODAY = date.today().isoformat()
DATE = config.get("date", TODAY)
