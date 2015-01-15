import settings
import os
import shutil

if __name__ == "__main__":
    path = settings.ETERNABOT_PATH + "strategies/"
    for f in os.listdir(path):
        if f[-3:] != ".py":
            continue
        current_file = path + f
        backup = current_file + ".bak"
        shutil.copy(current_file, backup)


