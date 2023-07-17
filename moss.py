import mosspy

userid = 350964312

m = mosspy.Moss(userid, "cc")

m.addFile("E:/YandexDisk/Archie/bio 2021/fiit/3/3/extra/2018_glazatov.txt")
m.addFile("E:/YandexDisk/Archie/bio 2021/fiit/3/3/extra/2020_shirokov.txt")
m.addFile("E:/YandexDisk/Archie/bio 2021/fiit/3/3/extra/2021_rudoy.txt")
m.addFile("E:/YandexDisk/Archie/bio 2021/fiit/3/3/extra/2021_tsarev.txt")

url = m.send(lambda file_path, display_name: print('*', end='', flush=True))
print()

print ("Report Url: " + url)

m.saveWebPage(url, "report.html")

mosspy.download_report(url, "report/", connections=8, log_level=10, on_read=lambda url: print('*', end='', flush=True))
