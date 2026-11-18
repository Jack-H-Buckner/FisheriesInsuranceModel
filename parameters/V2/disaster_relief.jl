
# probabilties of disaster responce
μ_TR = 3.8 # mean time until relief 
# this probabiltiy determines the mean time to relief 
# the relief is automatically delayed one year 
# (the event must occur before a payotu occurs)
# and then a geometric disturibution with mean 
# μ_TR-1 determines when the payout comes
pR = 1/(μ_TR-1) # probabiltiy of relief given decaration status 
pD = 0.9 # probablity a loss event is declared a disaster