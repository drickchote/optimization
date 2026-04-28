import csv
import matplotlib.pyplot as plt

FILE_NAME = "pareto_evolution.csv"

# Estrutura:
# data[generation] = [(risk, expectedReturn), ...]
data = {}

with open(FILE_NAME, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        gen = int(row["generation"])
        risk = float(row["risk"])
        variance = risk * risk
        expected_return = float(row["expectedReturn"])

        if gen not in data:
            data[gen] = []

        data[gen].append((variance, expected_return))

generations = sorted(data.keys())

if not generations:
    raise ValueError("Nenhuma geração encontrada no CSV.")

current_index = 0

fig, ax = plt.subplots()

def draw_generation():
    ax.clear()

    gen = generations[current_index]
    points = sorted(data[gen], key=lambda p: p[0])  # ordena por risk

    risks = [p[0] for p in points]
    returns = [p[1] for p in points]

    ax.plot(risks, returns, marker="o")
    ax.set_title(f"Pareto Frontier - Generation {gen}")
    ax.set_xlabel("Risk")
    ax.set_ylabel("Expected Return")
    ax.grid(True)

    fig.canvas.draw()

def on_key(event):
    global current_index

    if event.key == "right":
        if current_index < len(generations) - 1:
            current_index += 1
            draw_generation()

    elif event.key == "left":
        if current_index > 0:
            current_index -= 1
            draw_generation()

fig.canvas.mpl_connect("key_press_event", on_key)

draw_generation()
plt.show()