import sys
import time
import pandas as pd
import subprocess


def output_shell(line):
    """
    Save output of shell line that has pipes
    taken from: https://stackoverflow.com/questions/7389662/link-several-popen-commands-with-pipes
    :param line:
    :return:
    """
    try:
        shell_command = subprocess.Popen(line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    except OSError:
        return None
    except ValueError:
        return None

    (output, err) = shell_command.communicate()
    shell_command.wait()
    if shell_command.returncode != 0:
        print("Shell command failed to execute")
        return None

    return output

def du(path):
    return subprocess.check_output(['du', '-sh', path]).split()[0].decode('utf-8')

def bjobs():
    return output_shell("bjobs | grep -c 'RUN'")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else '.'

    bjobs_exist = True
    folder_size = []
    number_jobs = []
    total_time = []
    current_time = 0
    time_interval = 10

    try:
        with open('./monitoring.csv', 'a') as file:
            file.write('current_time,folder_size,num_bjobs\n')
            while True and bjobs_exist:
                folder_size.append(du(path))
                number_jobs.append(bjobs())
                total_time.append(current_time)
                to_write = str(current_time)+','+str(folder_size[-1])+','+str(number_jobs[-1])+'\n'
                file.write(to_write)
                current_time += time_interval
                time.sleep(time_interval)
                # if "No unfinished job found" in output_shell("bjobs"):
                #     bjobs_exist = False
    except KeyboardInterrupt:
        #time.sleep(time_interval)
        file.close()
        # d = {"folder_size": folder_size, "current_time": total_time, "num_bjobs": number_jobs}
        # df = pd.DataFrame(d)
        # df.to_csv("./monitoring.csv")
        raise