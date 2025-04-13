import pandas as pd

def calculate_centers_all_models(pdbqt_lines):
    results = []
    coords = []
    model_num = None
    in_model = False

    for line in pdbqt_lines:
        if line.startswith("MODEL"):
            model_num = int(line.split()[1])
            coords = []
            in_model = True
        elif line.startswith("ENDMDL"):
            if coords:
                xs, ys, zs = zip(*coords)
                center = (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
                results.append({
                    "Model": model_num,
                    "Center X": round(center[0], 3),
                    "Center Y": round(center[1], 3),
                    "Center Z": round(center[2], 3),
                })
            in_model = False
        elif in_model and line.startswith("HETATM"):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                continue  # Skip malformed lines

    return pd.DataFrame(results)

with open("/Users/alexbotkova/analysis_of_docking_studies/test_files/argatroban/result-2025-04-13T13_10_31.558Z/out_vina.pdbqt") as f:
    lines = f.readlines()

centers_df = calculate_centers_all_models(lines)
print(centers_df)
