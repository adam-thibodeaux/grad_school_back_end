from io import StringIO
import pandas as pd
def main(string, output_path):
    df = pd.read_csv(StringIO(string), sep="\s\s+")
    df.to_csv(output_path, sep=",")


xds_string = """RESOLUTION LIMIT   NUMBER OF REFLECTIONS OBSERVED    NUMBER OF REFLECTIONS UNIQUE     NUMBER OF REFLECTIONS POSSIBLE     COMPLETENESS OF DATA   R-FACTOR OBSERVED   R-FACTOR EXPECTED   R-FACTOR COMPARED   I/SIGMA       R-meas     CC(1/2)     Anomal      SigAno       Nano
     3.19         118      48       108       44.4%       9.3%     21.4%      118    2.25     12.2%    99.0*   -15    0.407       9
     2.26         219      78       172       45.3%      20.8%     25.3%      216    1.81     26.8%    97.8*   -13    0.750      28
     1.84         258      93       215       43.3%      26.6%     26.4%      254    1.99     33.4%    94.2*    23    0.852      33
     1.60         330     116       250       46.4%      42.4%     36.0%      321    1.70     52.9%    83.3*   -23    0.660      46
     1.43         383     136       286       47.6%      49.4%     50.2%      370    1.17     60.6%    86.7*     5    0.718      57
     1.30         377     132       295       44.7%      85.5%     78.8%      367    1.02    105.2%    54.4*     3    0.732      51
    total        1685     603      1326       45.5%      22.4%     27.9%     1646    1.53     28.3%    97.9*    -1    0.721     224"""
if __name__ == "__main__":
    main(xds_string, "./temp.csv")