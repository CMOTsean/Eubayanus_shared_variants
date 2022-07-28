import Filter_vcfs.filter_vcf as vcf
import pandas as pd
import matplotlib.pyplot as plt
import sys


chrom_dict = {"tig00000018_np1212121212_208603": 208603,
"tig00000020_np1212121212_1284184":        1284184,
"tig00000017_np1212121212_287396": 287396,
"tig00000014_np1212121212_424934": 424934,
"tig00000007_np1212121212_768921": 768921,
"tig00000011_np1212121212_766035": 766035,
"tig00000013_np1212121212_570983": 570983,
"tig00000008_np1212121212_837350": 837350,
"tig00000024_np1212121212_979429": 979429,
"tig00000003_np1212121212_1002340":        1002340,
"tig00000009_np1212121212_763993": 763993,
"tig00000015_np1212121212_319450": 319450,
"tig00000006_np1212121212_914339": 914339,
"tig00000002_np1212121212_1108862":        1108862,
"tig00000010_np1212121212_662407": 662407,
"tig00000004_np1212121212_1081690":        1081690,
"tig00000019_np1212121212_69556":  69556}

chrom_names = {"tig00000018_np1212121212_208603": 'Chr I',
"tig00000020_np1212121212_1284184":        'Chr II-IV',
"tig00000017_np1212121212_287396": 'Chr VI',
"tig00000014_np1212121212_424934": 'Chr IX',
"tig00000007_np1212121212_768921": 'Chr XIV',
"tig00000011_np1212121212_766035": 'Chr X',
"tig00000013_np1212121212_570983": 'Chr V',
"tig00000008_np1212121212_837350": 'Chr VII-XV',
"tig00000024_np1212121212_979429": 'Chr XIII',
"tig00000003_np1212121212_1002340":        'Chr IV-II',
"tig00000009_np1212121212_763993": 'Chr XV-VIII',
"tig00000015_np1212121212_319450": 'Chr III',
"tig00000006_np1212121212_914339": 'Chr XVI',
"tig00000002_np1212121212_1108862":        'Chr VII',
"tig00000010_np1212121212_662407": 'Chr XI',
"tig00000004_np1212121212_1081690":        'Chr XII',
"tig00000019_np1212121212_69556":  'mtDNA'}

def check_vars(s):
    if s[:3] in ["0/0", "0|0"]:
        return "homref"
    elif s[:3] in ["./.", ".|."]:
        return "nocall"
    else:
        return "var"

tibet_count=0
ucd_count=0
nc_count=0

fin = sys.argv[1]
out_dir = sys.argv[2]

fileout="blank.vcf"
df = pd.read_csv(fin, skiprows=vcf.skip_rows_check(fin, fileout), sep="\t")
df = df.rename(columns={"#CHROM": "CHROM"})

length_dict = {'UCD': 0,'tibet': 0,'NC': 0,'nocall': 0} #used later for calculating genome proportion


for chrom in chrom_dict:
    chrom_df = df[df["CHROM"] == chrom]
    past="W34-70"
    less_UCD=["yHRVM108", "yHRVM107", "new_CDFM21L", "ABFM5L.1"]
    less_tibet=["yHRVM108", "yHRVM107", "UCD646", "UCD650"]
    less_NC=["UCD646", "UCD650", "new_CDFM21L", "ABFM5L.1"]
    shared_vars=[]

    def find_shared_vars(row):
        """
        Takes a dataframe row as input, and retuurns whether the variant in 'past' is shared with Tibet, UCD or NC population (or is nocall)

        """
        global tibet_count
        global ucd_count
        global nc_count
        if check_vars(row[past]) == "nocall":
            return "nocall"
        elif (check_vars(row[past]) == "var") and (check_vars(row["new_CDFM21L"]) == "var"):
            for s in less_tibet:
                if check_vars(row[s]) == "var":
                    return "null"
            tibet_count+=1
            return "tibet"
        elif (check_vars(row[past]) == "var") and (check_vars(row["yHRVM107"]) == "var"):
            for s in less_NC:
                if check_vars(row[s]) == "var":
                    return "null"
            nc_count+=1
            return "NC"
        elif (check_vars(row[past]) == "homref") and (check_vars(row["UCD646"]) == "homref"):
            tot =0
            for s in less_UCD:
                if check_vars(row[s]) == "var":
                    tot += 1
            if tot == 4:
                ucd_count += 1
                return "UCD"
            else:
                return "null"
        else:
            return "null"
    
    shared_vars = chrom_df.apply(find_shared_vars, axis=1)
    plot_df = chrom_df[["POS"]].copy()
    plot_df["type"] = shared_vars
    plot_df = plot_df[plot_df["type"] != "null"]
    col_list = []
    start_list=[]
    stop_list = []
    for i in plot_df["type"]:
        if i == "nocall":
            col_list.append("k")
            start_list.append(0.13)
            stop_list.append(0.16)
        elif i == "UCD":
            col_list.append("#FF007F")
            start_list.append(0.07)
            stop_list.append(0.10)
        elif i == "tibet":
            col_list.append("#007FFF")
            start_list.append(0.04)
            stop_list.append(0.07)
        elif i == "NC":
            col_list.append("green")
            start_list.append(0.10)
            stop_list.append(0.13)
    
    #Create blocks from variants for plotting
    list_of_winds=[]
    cur_start = 0
    cur_pos = 0
    single_count=0
    cur_type = ""
    col_dict={"UCD": "#FF007F",
              "tibet": "#007FFF",
              "NC": "green",
              "nocall": "k"
              }
###    for l in plot_df[plot_df["type"] != "nocall"].itertuples():  #used to exclude nocalls from plot
    for l in plot_df.itertuples():
        if cur_type == "":
            cur_pos = l.POS
            cur_type = l.type
        elif l.type != cur_type and single_count == 0:
            list_of_winds.append([cur_start, round((cur_pos+l.POS)/2), col_dict[cur_type]])
            cur_start = round((cur_pos+l.POS)/2)
            cur_type = l.type
            cur_pos = l.POS
            single_count = 1
        elif l.type != cur_type and single_count == 1:
            cur_pos = l.POS
            cur_type = l.type
            single_count = 0
        else:
            cur_pos = l.POS
            single_count = 0
    list_of_winds.append([cur_start, chrom_dict[chrom],col_dict[cur_type]])
    windf = pd.DataFrame(list_of_winds, columns=["start", "end", "type"])
    for i,r in windf.iterrows():
        if r['type'] == 'k':
            length_dict['nocall'] += (r["end"]-r["start"])
        elif r['type'] == 'green':
            length_dict['NC'] += (r["end"] - r["start"])
        elif r['type'] == '#FF007F':
            length_dict['UCD'] += (r["end"] - r["start"])
        elif r['type'] == '#007FFF':
            length_dict['tibet'] += (r["end"] - r["start"])

    
    #Set up figure
    plt.rc('font', size=38)
    fig, ax = plt.subplots()
    fig.set_size_inches(50, 5)

    #Set up plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    plt.ylim(-0.06,0.16)
    plt.xlim(-20, chrom_dict[chrom])
    plt.ticklabel_format(style='plain')
    plt.yticks([])
    plt.xlabel("Position (bp)")
    plt.subplots_adjust(bottom=0.25)
    
    #Plot and save data
    plt.vlines(ymin=start_list, ymax=stop_list, x=plot_df["POS"], colors=col_list, linewidth=3)
    for i in windf.itertuples():
        plt.hlines(xmin=i.start, xmax=i.end, y=0, colors=i.type, linewidth=60)
    fig.savefig(out_dir+'/'+past+"_"+chrom_names[chrom]+".png", dpi=400)
    plt.close()

#print("Tibet Count", str(tibet_count))
#print("UCD Count", str(ucd_count))
#print("NC Count", str(nc_count))
