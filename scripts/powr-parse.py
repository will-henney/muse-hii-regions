from pathlib import Path

datapath = Path("../stars/powr")
datafiles = datapath.glob("*_colors.txt")

rslt = [
    ["Grid", "Model", "Q0", "Q1", "Q2", "Q3"],
    None,
]

for datafile in datafiles:
    with open(datafile) as f:
        lines = f.readlines()

    _, grid_id, _, model_id = lines[0].split()
    # print(grid_id, model_id)
    _, _, _, Q0 = lines[6].split()
    _, _, _, Q1 = lines[7].split()
    _, _, _, Q2 = lines[8].split()
    _, _, _, Q3 = lines[9].split()

    rslt.append([grid_id, model_id, Q0, Q1, Q2, Q3])

print(rslt)
