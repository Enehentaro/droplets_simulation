import os
import pandas as pd
import glob
import matplotlib.pyplot as plt


def patientAverage(df:pd.DataFrame) -> pd.DataFrame:
    df_s = df.sort_index()
    # print(df_s)

    num_operation:int = len(df_s) // 5

    df_new = None

    for i in range(num_operation):
        df_inCase = df_s.iloc[5*i : 5*(i+1)]
        s_mean = df_inCase.mean()

        casename:str = df_s.index[5*i]
        casename = casename.rsplit("_", 1)[0]

        s_mean["casename"] = casename
        df_mean = pd.DataFrame([s_mean])
        if df_new is None:
            df_new = df_mean
        else:
            df_new = pd.concat([df_new, df_mean])

    df_new.set_index("casename", inplace=True)
    # print(df_new)
    return df_new

    # while len(index_list) > 0:
    #     casename:str = index_list[0]
    #     casename = casename.rsplit("_", 1)[0]
    #     case_list = [casename + "_" + patient for patient in ["A", "B", "C", "D", "E"]]
    #     df_inCase = df.loc[case_list]

    #     if df_new is None:
    #         df_new = df_inCase
    #     else:
    #         df_new = pd.concat([df_new, df_inCase])

    #     print(df_new)
    #     for case_remove in case_list:
    #         index_list.remove(case_remove)

def plot_var():
    for startSecond in startSecond_list:
        fname = f"./CountResults/sitting/count_from{startSecond}sec.csv"
        df_sitting = pd.read_csv(fname)

        fname = f"./CountResults/standing/count_from{startSecond}sec.csv"
        df_standing = pd.read_csv(fname)

        df_total = df_sitting #+ df_standing

        df_total["RoI_norm"] = df_total["RoI"] / df_total["RoI"].max()

        total_var_list.append(df_total["RoI_norm"].var())

        #####################

        fname = f"./CountResults/sitting/count_from{startSecond}sec_patientAverage.csv"
        df_sitting = pd.read_csv(fname)

        fname = f"./CountResults/standing/count_from{startSecond}sec_patientAverage.csv"
        df_standing = pd.read_csv(fname)

        df_patientAverage = df_sitting #+ df_standing

        df_patientAverage["RoI_norm"] = df_patientAverage["RoI"] / df_patientAverage["RoI"].max()

        patientAverage_var_list.append(df_patientAverage["RoI_norm"].var())


if __name__ == "__main__":
    startSecond_list = [i+1 for i in range(10)]

    total_var_list = []
    patientAverage_var_list = []

    df_officeSize = pd.read_csv("OfficeSize.csv", index_col="officename")

    for startSecond in startSecond_list:

        # file_list = glob.glob(f"count_results_sitting_221018/*from{startSecond}sec*")
        file_list = glob.glob(f"count_results_standing_221018/*from{startSecond}sec*")
        print(file_list)

        df_total = None
        df_patientAverage = None

        for filename in file_list:
            bname = os.path.basename(filename)
            officename = bname.split("_")[0]
            df_office = df_officeSize.loc[officename]
            office_area = df_office["Lx"] * df_office["Ly"]

            df_read = pd.read_csv(filename, index_col="casename")
            df_read.rename(index=lambda s: officename + "_" + s, inplace=True)
            df_read["RoI"] = df_read["num_drop"] / office_area

            if df_total is None:
                df_total = df_read
            else:
                df_total = pd.concat([df_total, df_read])

            df_ave = patientAverage(df_read)
            # print(df_ave)
            if df_patientAverage is None:
                df_patientAverage = df_ave
            else:
                df_patientAverage = pd.concat([df_patientAverage, df_ave])

        # outputDir = "countResults_sitting"
        outputDir = "countResults_standing"
        os.makedirs(outputDir, exist_ok=True)
        df_total.to_csv(outputDir + f"/count_from{startSecond}sec.csv")
        df_patientAverage.to_csv(outputDir + f"/count_from{startSecond}sec_patientAverage.csv")

        total_var_list.append(df_total["RoI"].var())
        patientAverage_var_list.append(df_patientAverage["RoI"].var())

    # print(var_list)


    #グラフプロット準備
    plt.plot(startSecond_list, total_var_list, label="total")
    plt.plot(startSecond_list, patientAverage_var_list, label="patient average")

    plt.title("Relation between StartTime and Variance")
    plt.xlabel("start time [sec]")
    plt.ylabel("variance")

    plt.xlim(0,10)

    plt.legend() #凡例を表示

    #グラフ出力
    plt.show()