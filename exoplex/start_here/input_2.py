
n_mod    = 2           # how many models will you run?
indp     = 'm'         # set mass or radius as independent value. M=mass, R=radius (default = M)
X        = [1,1.05]        # list of value(s) of independent parameter (M in ME=5.972e24 kg and R in RE=6.371e6 m 
fix_core = False       # fix core mass fraction? (True or False)    
cmf      = [0.32,0.]     # list core mass fraction for each model if fix_core = True (0,1)
Si_wt    = [0.0,0.12]       # wt fraction of core made of Si (0,<1) 
O_wt     = [0.0,0.0]       # wt fraction of core made of O  (0,<1)
S_wt     = [0.0,0]       # wt fraction of core made of Si (0,<1)
FeMg     = [1.,1.2]       # list elemental Fe/Mg ratios for whole planet or just mantle if fix_core = True 
SiMg     = [1.,1.3]       # list elemental Si/Mg ratios for whole planet or JUST THE MANTLE if fix_core = True                     
CaMg     = [0.3,.2]       # list elemental Ca/Mg ratios for whole planet or JUST THE MANTLE if fix_core = True     
AlMg     = [.03,.3]       # list elemental Al/Mg ratios for whole planet or JUST THE MANTLE if fix_core = True
xFeO     = [0.1,.3]      # list, Fraction of total Fe (mantle+core) between **0-1** in the mantle  [default = [0,] if fix_core = True]
wt_frac_water = [0.0,0.03] # list water mass fraction
#--------------------------------------------------------------------------------------#
