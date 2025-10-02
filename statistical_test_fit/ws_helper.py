def set_constant(w):
    print("--> Setting all var as constants")
    for v in w.allVars().contentsString().split(","):
        if v != "boson_mass" and v != "upsilon_mass" and v != "evt_weight":
            print("{v} to constant...".format(v=v))
            w.var(v).setConstant()
            w.var(v).Print()
    return w
