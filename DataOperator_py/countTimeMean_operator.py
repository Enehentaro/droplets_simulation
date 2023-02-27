import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import sys

"""
時系列データを操作し、カウントの時間平均値を算出してファイル出力するスクリプト
"""

def transpose(file_list:list, output_dir:str):
    """CSVファイルの行列を転置して、別ファイルに保存

    Args:
        file_list (list): CSVファイルリスト
        output_dir (str): 出力ディレクトリ
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in file_list:
        df_read = pd.read_csv(file_name, index_col=0, header=None)

        bname = os.path.basename(file_name)

        df_read.T.to_csv(output_dir+bname, index=False)
    # print()
    

def patientMean(file_list:list, output_dir:str):
    """感染者5人の平均を求め、別ファイルに保存

    Args:
        file_list (list): ファイルリスト
        output_dir (str): 出力ディレクトリ
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    patient_list = ["A", "B", "C", "D", "E"]

    for file_name in file_list:
        df_read = pd.read_csv(file_name)
        column_list = list(df_read.columns)

        bname = os.path.basename(file_name)

        df_new = pd.DataFrame()

        #このリストが空になるまでループ
        while len(column_list) > 0:

            column:str = column_list[0]
            case_name = column.rsplit("_", 1)[0]
            
            case_list = [case_name+"_"+patient for patient in patient_list]
            
            df_case = df_read[case_list]    #感染者5人分の列だけ抽出
            # print(df_case)
            s_mean = df_case.mean(axis=1)   #感染者5人の平均をとる
            # print(s_mean)

            df_new[case_name] = s_mean     #平均値をメインのDataFrameに追加

            column_list = [i for i in column_list if i not in case_list]    #今処理した5ケースをリストから削除

        df_new.sort_index(axis=1, inplace=True)
        # print(df_new)
        df_new.to_csv(output_dir+bname, index=False)
    # print()

def make_summaryCSV(file_list:list, output_fname:str):
    """別々に保存されている各オフィスのデータをひとつのファイルにまとめる

    Args:
        file_list (list): ファイルリスト
        output_fname (str): 出力ファイル名
    """

    df_main = pd.DataFrame()

    for file_name in file_list:
        df_read = pd.read_csv(file_name)

        bname = os.path.basename(file_name)
        officename = bname.split(".")[0]
        df = df_read.rename(columns=lambda s: officename + "_" + s)

        print(df)
        # print(f"{bname}: {df.shape}")

        df_main = pd.concat([df_main, df], axis=1)

    df_main.sort_index(axis=1, inplace=True)
    print(df_main)
    df_main.to_csv(output_fname, index=False)

def subtractDataFrame(macro_file_list:list, output_dir:str):
    """全体観測データから、机の上観測データを引き算して、新規ファイルに保存

    Args:
        macro_file_list (list): 全体観測結果ファイルリスト
        output_dir (str): 出力ディレクトリ
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for macro_file_name in macro_file_list:
        macro_file_name:str
        df_macro = pd.read_csv(macro_file_name)

        aboveDesk_file_name = macro_file_name.replace("fromSittingToStanding", "fromSittingToStanding_aboveDesk")
        df_aboveDesk = pd.read_csv(aboveDesk_file_name)

        df_exceptDesk = df_macro - df_aboveDesk

        bname = os.path.basename(macro_file_name)

        negative_count = (df_exceptDesk < 0).values.sum()
        if negative_count >= 1:
            # 負の値が混じっていたらエラー
            sys.stderr.write(f"NegativeValueFound: {negative_count} in `{bname}`\n")

        df_exceptDesk.to_csv(output_dir+bname, index=False)

if __name__ == "__main__":
    base_dir = "Count_timeSeries/onlyFloating/"
    # base_dir = "Count_timeSeries/fromSittingToStanding_aboveDesk/"

    file_list = glob.glob(base_dir+"office*.csv")
    transpose(file_list, output_dir=base_dir+"transpose/")

    # file_list = glob.glob("Count_timeSeries/fromSittingToStanding/transpose/office*.csv")
    # subtractDataFrame(file_list, output_dir=base_dir+"transpose/")

    file_list = glob.glob(base_dir+"transpose/office*.csv")
    patientMean(file_list, output_dir=base_dir+"patientMean/")

    file_list = glob.glob(base_dir+"patientMean/office*.csv")
    make_summaryCSV(file_list, output_fname=base_dir+"patientMean/summary.csv")


    # 時系列データを可視化
    # # df = pd.read_csv("Count_timeSeries/exceptDesk/patientMean/summary.csv")
    # df = pd.read_csv("Count_timeSeries/onlyFloating/transpose/office3.csv")
    # df["960_246_aout_E"].plot(ylim=(0,10000))
    # plt.xlabel("Time Step")
    # plt.ylabel("Count")
    # # df.iloc[1:, :20].plot(ylim=(0,10000), figsize=(16,9))
    # # df_fil = df.filter(like="office1_", axis=1)
    # # df_fil.iloc[1:].plot(ylim=(0,10000), figsize=(16,9))
    # plt.show()

    # 時系列データを平均化（時間平均）
    # s_mean = df.mean()
    # print(s_mean)
    # print(s_mean.var())



# ==========================================================

# def patientAverage(ser:pd.Series) -> pd.Series:

#     num_operation:int = len(ser) // 5

#     s_new = None

#     for i in range(num_operation):
#         s_inCase = ser.iloc[5*i : 5*(i+1)]
#         mean = s_inCase.mean()

#         casename:str = ser.index[5*i]
#         casename = casename.rsplit("_", 1)[0]

#         s_mean = pd.Series({casename : mean})
#         if s_new is None:
#             s_new = s_mean
#         else:
#             s_new = pd.concat([s_new, s_mean])

#     # print(df_new)
#     return s_new.sort_index()


# ==========================================================


# import pandas as pd
# import matplotlib.pyplot as plt


# df_read = pd.read_csv("Count_timeSerirs/sitting/office1.csv", header=None)

# df = df_read.rename(columns={0: 'case_name'})

# df.set_index('case_name', inplace=True)

# df = df.T

# print(df)

# # df.iloc[:,:10].plot()
# df[["0_0_A", "0_0_B", "0_0_C", "0_0_D", "0_0_E"]].plot(figsize=(9, 6))


# # これがないとビューアーが起動しない（VSCodeだけ？）
# plt.show()





# ==========================================================
# import glob

# file_list = glob.glob("Count_timeSerirs/standing/*")
# print(file_list)

# for file_name in file_list:

#     with open(file_name, 'r') as f:
#         text_list = f.readlines()

#     new_file_name = file_name.replace("CountTimeSeries_standing_office", "office")
#     new_file_name = new_file_name.replace("sec.csv", ".csv")
#     with open(new_file_name, "w") as f:
#         f.write("".join(text_list[1:]))

# ==========================================================

# import glob
# import os

# file_list = glob.glob("Count_timeSerirs/sitting/*")
# print(file_list)

# for file_name in file_list:

#     new_file_name = file_name.replace("sec.csv", ".csv")
#     os.rename(file_name, new_file_name)
