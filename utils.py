def pretty_print_prediction(prediction):
  for aa, state in prediction:
    s = ""
    if state == 0:
      s = "Variable Region"
    if state in range(1,11):
      s = "Beta Strand"
    if state == 11:
      s = "Variable Region"
    if state in [12,13]:
      s = "Binding Residues"
    if state == 14:
      s = "Variable Region"
    if state in range(15,24):
      s = "Conserved SDPYVK Region"
    if state == 24:
      s = "Variable Region"
    if state in range(25,36):
      s = "Beta Strands"
    if state == 36:
      s = "Variable Region"
    if state in range(37,58):
      s = "Final Beta Strand in model" 
    if state == 58:
      s = "Variable Region"
    print s + " " + aa
