# Importar librerías
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Funciones
def complement(genome):
    complementary= {"A":"T", "C":"G", "G":"C", "T":"A"}
    c_genome =""
    for letter in genome:
        c_genome += complementary[letter]
    return c_genome

def reverse(genome):
    return genome[::-1]

def reverse_complement(genome):
    c_genome = complement(genome)
    reverse_complement_genome= reverse(c_genome)
    return reverse_complement_genome

def similarities(genome, gene, step_size=1):
    window_size = len(gene)
    iter_limit = len(genome) - window_size + 1
    positions = []
    similarities = []
    for x in range(0, iter_limit, step_size):
        window = genome[x:x + window_size]
        positions.append(x)
        matches = sum(1 for i in range(len(window)) if window[i] == gene[i])
        similarity = matches / window_size
        similarities.append(similarity)
    return positions, similarities


# Variables iniciales
genome_file = "solanum_phureja_chromosome1.fasta"
cds_file = "CDS_ribosomal_protein_L2_cromosome1_S_phureja.fasta"
genome = ""
gene = ""
gene_name = ""
genome_name = ""

# Leer archivos
with open(genome_file) as file:
    for line in file:
        if line.startswith(">"):
            genome_name = line.strip()
        else:
            genome += line.strip()
with open(cds_file) as file2:
   for line in file2:
        if line.startswith(">"):
            gene_name = line.strip()
        else:
            gene += line.strip()

# Calcular el  reverso complementario del gen
gene_reverse = reverse_complement(gene)

# Comparacion del CDS con el genoma
positions_forward, similarities_forward = similarities(genome, gene)
positions_reverse, similarities_reverse = similarities(genome, gene_reverse)

# Crear DataFrame
df = pd.DataFrame({"Position": positions_forward,
                     "Similitud_Forward": similarities_forward,
                     "Similitud_Reverse": similarities_reverse})
genome_name = "_".join(genome_name.split()[1:3])
df_name = f"similitud_with_{genome_name}.csv"
df.to_csv(df_name)

# # Importar DataFrame
# df = pd.read_csv("similitud_with_Solanum_phureja.csv")
# position_real = df["Position"] + 1
# df["Position"] = position_real
# df = df[["Position", "Similitud_Forward", "Similitud_Reverse"]]

# Convertir DataFrame de formato ancho a largo
data = pd.melt(df, id_vars="Position", value_vars=["Similitud_Forward", "Similitud_Reverse"], var_name= "Chain", value_name ="Similitud")


# Grafica 1
sns.lineplot(data=df, x="Position", y="Similitud_Forward")
plt.grid()
plt.xlabel("Posición en el Genoma")
plt.ylabel("Similitud Forward")
plt.title("Similitud Forward entre el CDS y el Genoma")
plt.savefig("1_Similitud_CDS_Forward_Genoma.pdf")
plt.show()

# Grafica 2
sns.lineplot(data=df, x="Position", y="Similitud_Reverse")
plt.grid()
plt.xlabel("Posición en el Genoma")
plt.ylabel("Similitud Reverse")
plt.title("Similitud Reverse entre el CDS y el Genoma")
plt.savefig("2_Similitud_CDS_Reverse_Genoma.pdf")
plt.show()

# Grafica 3
sns.lineplot(data=df, x="Position", y="Similitud_Forward", label="Similitud Forward", color="blue")
sns.lineplot(data=df, x="Position", y="Similitud_Reverse", label="Similitud Reverse", color="purple")
plt.xlabel("Posición en el Genoma")
plt.ylabel("Similitud")
plt.title("Comparación de Similitudes Forward y Reverse")
plt.grid()
plt.savefig("3_Similitud_CDS_Forward_Reverse_Genoma.pdf")
plt.show()

# Grafica 4
max_forward = df["Similitud_Forward"].max()
max_reverse = df["Similitud_Reverse"].max()

sns.lineplot(data=df, x="Position", y="Similitud_Forward", label="Similitud Forward", color="blue")
sns.lineplot(data=df, x="Position", y="Similitud_Reverse", label="Similitud Reverse", color="red")
plt.axhline(max_forward, color="blue", linestyle="--", label=f"Máximo Forward: {max_forward}")
plt.axhline(max_reverse, color="red", linestyle="--", label=f"Máximo Reverse: {max_reverse}")
plt.xlabel("Posición")
plt.ylabel("Similitud")
plt.title("Similitudes con Máximos Resaltados")
plt.legend()
plt.grid()
plt.savefig("4_Similitud_CDS_Forward_Reverse_Genoma_Max.pdf")
plt.show()

# Grafica 5
threshold = 0.30
highlight_positions = df[df["Similitud_Forward"] > threshold]

sns.lineplot(data=df, x="Position", y="Similitud_Forward", label="Similitud Forward", color="blue")
plt.scatter(highlight_positions["Position"], highlight_positions["Similitud_Forward"], color="orange", label="Puntos Destacados")
plt.xlabel("Posición")
plt.ylabel("Similitud Forward")
plt.title(f"Similitud Forward con Umbral {threshold}")
plt.legend()
plt.grid()
plt.savefig("5_Similitud_CDS_Forward_Reverse_Genoma_threshold.pdf")
plt.show()

# Grafica 6
sns.set_theme(style="whitegrid")
sns.histplot(data= data,
            x="Similitud",
            hue="Chain",
            palette="dark",
            alpha=0.6,)
plt.xlim(0.2, 0.3)
plt.xlabel("Similitud")
plt.ylabel("Frecuencia")
plt.title("Distribución Comparativa de Similitudes Forward y Reverse")
plt.savefig("6_Distribucion_Similitudes_Forward_Reverse.pdf")
plt.show()

# Grafica 7
heatmap_data = df[["Similitud_Forward", "Similitud_Reverse"]].T
custom_palette = sns.dark_palette("#69d", reverse=True, as_cmap=True)
sns.heatmap(heatmap_data,cmap= custom_palette, vmin= 0.2, vmax = 0.9, cbar=True)
plt.title("Mapa de Calor de Similitudes")
plt.xticks(rotation = 30)
plt.xlabel("Posiciones en el Genoma")
plt.ylabel("Similitud")
plt.savefig("7_MapaCalor_Similitudes_Forward_Reverse.pdf")
plt.show()

#  Grafica 8
heatmap_data = df[["Similitud_Forward"]].T
custom_palette = sns.dark_palette("#69d", reverse=True, as_cmap=True)
sns.heatmap(heatmap_data,cmap= custom_palette, vmin= 0.2, vmax = 0.9, cbar=True)
plt.title("Mapa de Calor de Similitudes Forward")
plt.xticks(rotation = 30)
plt.xlim(37000,40000)
plt.xlabel("Posiciones en el Genoma")
plt.ylabel("Similitud")
plt.savefig("8_MapaCalor_Similitudes_Forward.pdf")
plt.show()

# Grafica 9
sns.set(style="whitegrid")
sns.scatterplot(data=df, x="Position", y="Similitud_Forward", color="purple")
plt.xlabel("Posicion")
plt.ylabel("Similitud Forward")
plt.title("Similitud Forward vs Posicion")
plt.savefig("9_Similitud_scatterplot_Forward_Genoma.pdf")
plt.show()

# Grafica 10
sns.set(style="whitegrid")
sns.scatterplot(data=df, x="Position", y="Similitud_Reverse", color="purple")
plt.xlabel("Posicion")
plt.ylabel("Similitud Reverse")
plt.title("Similitud Reverse vs Posicion")
plt.savefig("10_Similitud_scatterplot_Reverse_Genoma.pdf")
plt.show()

# Grafica 11
sns.set(style="whitegrid")
sns.scatterplot(data=data,
                x="Position",
                y="Similitud",
                hue = "Chain")
plt.legend()
plt.xlabel("Posicion")
plt.ylabel("Similitud")
plt.title("Similitud Forward y Reverse vs Posicion")
plt.savefig("11_Similitud_scatterplot_Forward_Reverse_Genoma.pdf")
plt.show()

print("Fin")