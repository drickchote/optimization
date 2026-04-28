array = [0,0.0554943,1.89077,0,0.120709,0,0,0.0624356,0,0,0,0,0.166856,0,0.165253,0,0,0,0,0.0432313,0.0853781,0,0,0.0516612,0,0,0,0.0489074,0.00633353,0.104663,0]
normalized_array = []

lower_bound = 0.01
picked_count = 0

for i in array:
    if i > 0:
        picked_count+=1

def normalize(x):

    return lower_bound + x/sum(array) * (1 - lower_bound * picked_count)


for i in array:
    if i == 0:
        continue
    normalized_array.append(normalize(i))
    

print(sum(array))
print(sum(normalized_array))