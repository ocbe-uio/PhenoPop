import matplotlib
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import colorsys

font = {'family' : 'DejaVu Sans',
        'weight' : 'semibold',
        'size'   : 16}

matplotlib.rc('font', **font)

def desaturate(color, degree):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    c = colorsys.hls_to_rgb(c[0], 1 - degree * (1 - c[1]), c[2])
    return c

old_mixturecolors = ['tab:blue','r'] #, 'g']

magikarp = [
'#ffde29',
'#bd2052',
'#bd4141',
'#ff9c62',
'#830041',
'#7b6252',
'#c5ac73',
'#cdcdd5',
'#525262',
'#8b8ba4'
]

gyarados = [
'#41b4ee',
'#cdb47b',
'#186294',
'#184173',
'#737b7b',
'#5a4120',
'#6a1820',
'#ee6241',
'#bd3162'
]

starmie = [
'#8b73bd',
'#313973',
'#c59420',
'#5a4131',
'#eec54a',
'#cd205a',
'#8b1052',
'#eeb4cd'
]

together = [
'#313973',
'#41b4ee',
'#bd2052',
'#8b73bd',
'#ffde29'
]

homemade_8 = [
"#0094C6",
"#D64933",
"#F6AE2D",
"#27187E",
"#0C7C59"
]
homemade_9 = [
"#0094c6",
"#ff0035",
"#e89005",
]

# https://carto.com/carto-colors/
cartocolors_antique = [
"#855C75",
"#D9AF6B",
"#AF6458",
"#736F4C",
"#526A83",
"#625377",
"#68855C",
"#9C9C5E",
"#A06177",
"#8C785D",
"#467378",
"#7C7C7C"
]
cartocolors_bold = [
"#7F3C8D", # purple
"#3969AC", # blue
"#11A579", # green
"#008695", # teal blue
"#80BA5A", # light green
"#F2B701", # yellow 
"#E68310", # orange
"#E73F74", # pink
"#CF1C90", # sharp pink
"#f97b72", # pale pink
"#4b4b8f", # marine blue
"#A5AA99", # slate
]
cartocolors_pastel = [
"#66C5CC",
"#F6CF71",
"#F89C74",
"#DCB0F2",
"#87C55F",
"#9EB9F3",
"#FE88B1",
"#C9DB74",
"#8BE0A4",
"#B497E7",
"#D3B484",
"#B3B3B3"
]
cartocolors_prism = [
"#5F4690",
"#1D6996",
"#38A6A5",
"#0F8554",
"#73AF48",
"#EDAD08",
"#E17C05",
"#CC503E",
"#94346E",
"#6F4070",
"#994E95",
"#666666"
]
cartocolors_safe = [
"#88CCEE",
"#CC6677",
"#DDCC77",
"#117733",
"#332288",
"#AA4499",
"#44AA99",
"#999933",
"#882255",
"#661100",
"#6699CC",
"#888888"
]
cartocolors_vivid = [
"#E58606",
"#5D69B1",
"#52BCA3",
"#99C945",
"#CC61B0",
"#24796C",
"#DAA51B",
"#2F8AC4",
"#764E9F",
"#ED645A",
"#CC3A8E",
"#A5AA99"
]

cartocolors_3_5 =  ["#7f3c8d","#3969ac","#80ba5a","#cf1c90","#e73f74","#f97b72","#e68310","#f2b701"]
paul_tol =         ['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499', '#DDDDDD']
paul_tol =         ['#332288', '#882255', '#CC6677', '#117733', '#88CCEE', '#DDCC77', '#44AA99', '#999933', '#AA4499', '#DDDDDD'] # reordered
paul_tol_darker =  ["#2a1c71", "#711c47", "#c04156", "#0e652b", "#399082", "#52b5e6", "#d2bc4d", "#81812b", "#903982", "#bababa"] # heavier
paul_tol_lighter = ["#412bad", "#ad2b6c", "#d37b8a", "#18a647", "#59bdac", "#97d3f0", "#e2d48b", "#b9b93e", "#bd59ac", "#e2e2e2"] # lighter
# Changed to get some pink
paul_tol =         ["#332288", "#882255", "#cc6677", "#aa4499", "#44aa99", "#88ccee", "#ddcc77", "#999933", "#117733", "#dddddd"] # reordered again again
paul_tol_darker =  ["#261965", "#65193f", "#ac394d", "#813374", "#338174", "#35a8e2", "#ccb334", "#737326", "#0d5a27", "#a6a6a6"] # darker, brightness -25
paul_tol_lighter = ["#5037cd", "#cd3782", "#da8f9c", "#c772b9", "#72c7b9", "#a7d9f2", "#e6da9d", "#c7c757", "#1ed059", "#e7e7e7"] # lighter, brightness +27
limitationcolors = ['#117733', '#999933', '#AA4499', '#332288', '#DDCC77', '#88CCEE', '#44AA99', '#DDDDDD']
limitationcolors = ['#AA4499', '#999933', '#117733', '#332288', '#332288', '#88CCEE', '#44AA99', '#DDDDDD']
limitationcolors = ['#AA4499', '#999933', '#117733', '#332288', '#35a8e2', '#332288', '#44AA99', '#DDDDDD']
paul_tol =         ["#332288", "#882255", "#cc6677", "#aa4499", "#44aa99", "#88ccee", "#ddcc77", "#999933", "#117733", "#dddddd"] 

paul_tol =         ["#332288", "#cc6677", "#ddcc77"] # these are wrong:, "#aa4499", "#44aa99", "#88ccee", "#ddcc77", "#999933", "#882255", "#dddddd"] # reordered again**3
limitationcolors_case_2 = ["#332288", "#cc6677", "#117733"]
paul_tol_darker =  ["#261965", "#65193f", "#ac394d", "#813374", "#c04156", "#35a8e2", "#338174", "#737326", "#0d5a27", "#a6a6a6"] # darker, brightness -25
paul_tol_lighter = ["#5037cd", "#cd3782", "#da8f9c", "#c772b9", "#d37b8a", "#a7d9f2", "#72c7b9", "#c7c757", "#1ed059", "#e7e7e7"] # lighter, brightness +27

# For circles in elbow plots: 
# 332288f4
# 44aa99f4
# aa4499f4

greycolor = '#DDDDDD'

homemade_10 = ["#2364AA","#ff0035","#ffc71f","#00b8f5","#679436"]
homemade_12 = ["#d92622","#387fb6","#4caf49","#934e9a","#ec7c1c"]
naturepaper = ["#387fb6","#d92622","#4caf49","#934e9a","#ec7c1c","#d23a8d","#6abea4","#f08c65","#8984be","#dd9b2b"]
banks = ["#63b9ce","#e68447","#e6d257","#5b7e21","#ac5083","#6abea4","#f2e30b","#8984be", "#334494"]
colorbrewer2 = ['#7570b3','#d95f02','#1b9e77','#e7298a','#66a61e','#e6ab02','#a6761d']
seacolors = [together, homemade_8] #  gyarados, starmie]

def plot_moodboard(colorscheme):
    x = np.linspace(0,4*np.pi,1000)
    fig = plt.figure()
    for idx, color in enumerate(colorscheme):
        #plt.plot([0,1], [idx, idx], linewidth=20, color=color)
        plt.plot(x, np.sin(x + idx*np.pi/4), linewidth=2.5, color=color)
    #plt.title(schemename)
    
    #plt.set_axis_off()
    #plt.spines['top'].set_visible(False)
    #plt.spines['right'].set_visible(False)
    #plt.spines['bottom'].set_visible(False)
    #plt.spines['left'].set_visible(False)

    plt.show()

if __name__=="__main__":
    for colorscheme in seacolors:
        plot_moodboard(colorscheme)
