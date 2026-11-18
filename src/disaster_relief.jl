
# Update the disaster relief indicator. `ID` is the current declaration status (0/1),
# `IM_next` indicates a high mortality event this year (0/1), and IRrng is a draw that
# determines whether relief arrives. A disaster is always declared after a high
# mortality event (probability 1, so there is no random draw for it). If a disaster is
# already declared, relief `p.Relief` is paid with probability p.pR (when
# `p.pR > IRrng`) and the declaration ends; a new high mortality event then declares a
# disaster again, including in the year relief is paid. Returns the next declaration
# status and the relief payment.
function relief(ID,IM_next,IRrng,p)
    Rel = 0.0
    if (ID == 1) && (p.pR > IRrng)
        Rel = p.Relief
        ID = 0
    end 
    if (IM_next == 1) 
        ID = 1
    end 
    return ID, Rel
end 
