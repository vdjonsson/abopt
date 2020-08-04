import pandas as pd 
import seaborn as sb 
import matplotlib.pyplot as plt 
import antibody_pipeline as ab 


context = 'paper' 
style='ticks'
def plot_ddg_distribution(p1, p2, title):

    ddg = ab.calculate_ddg_bind(p1 , p2, title)
    sb.set(context=context, style='ticks')
    sb.distplot( ddg.ddg.values, kde=False, hist=True)
    plt.semilogy()
    plt.ylabel('log10(ddg)')
    plt.title(title)
    plt.show()


def plot_ddg_heatmap_structures(p1, p2, title):

    ddg = ab.calculate_ddg_bind(p1 , p2, title)

    ddg = ddg[['mut','ddg']]
    ddg = ddg.set_index('mut')
    ddg = ddg.transpose()

    sb.set(context=context, style=style)
    sb.heatmap(ddg, square=False)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_ddg_heatmap(ddg, title):

    ddg = ddg.transpose()

    sb.set(context=context, style=style)
    sb.heatmap(ddg, square=False)
    plt.title(title)
    plt.tight_layout()
    plt.show()

def combine_data_sets(ddg_array):

    merged = pd.DataFrame()
    merged = ddg_array[0]

    for ddg_data in ddg_array[1:]: 
        merged = merged.merge(ddg_data, on='mut') 
    
    merged = merged.drop_duplicates(subset=['mut'], keep='first')
    merged = merged.set_index('mut')
    return merged 

# mutations of RBD in context of C105_wt

# mutations of RBD in context of C105_mt  
p1 = '6XCM_Repair_TH28D_Repair'
p2 = '6XCM_Repair_TH28D_Repair_less_ab_Repair' #less ab, mutate virus 

#plot_ddg_distribution(p1, p2, title='ddg(mut(RBD)/6XCM_TH28D)')
#plot_ddg_heatmap(p1, p2, title='ddg(mut(RBD)/6XCM_TH28D)')


ddg1 = ab.calculate_ddg_bind(p1 , p2, ddg_name='ab1') 
ddg2 = ab.calculate_ddg_bind(p1 , p2, ddg_name='ab2') 
ddg3 = ab.calculate_ddg_bind(p1 , p2, ddg_name='ab3')
 
ddg_array = [ddg1, ddg2, ddg3]

ddgs = combine_data_sets(ddg_array)

plot_ddg_heatmap(ddgs, title='ddg(mut(RBD)/6XCM_TH28D)')



# postive values are destabilizing and therefore of concern 

