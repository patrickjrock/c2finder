from fasta import Fasta

def get_records(fname):
  f = open(fname,'r')
  records = []
  record = ""
  for line in f:
    if line[0] is '>' and record != "":
      records.append(record)
      record = "" 
    if line[0] != '>':
      record += line[:-1]
  return records

def count_emissions(labels, data):
  e = [[] for i in range(0,max(labels)]
  for i in range(0, len(labels)):
    if data[i] != '-':
      state = int(labels[i])
      print state
      e[state-1].append(data[i])
  return e

def emissions(fname, labels):
  """ emissions are sorted into bins for each state. The emissions for state 1 are in 
      emissions(fname)[0]"""
  records = get_records('fasta/out.fasta')
  es = []
  for i in range(0,len(records)):
    es.append(count_emissions(labels, records[i]))
  final = []
  for e in es:
    for i, state in enumerate(e):
      if i+1 > len(final):
        final.append([])
      final[i].extend(state)
  return final


def count_transitions(labels, data):
  t = []
  last = -1 # last observed state, used to do transition
  for i in range(0, len(labels)):
    if data[i] != '-':
      state = int(labels[i])
      if last != -1:
        t.append((last, state))
      last = state
  return t

def transitions(fname, labels):
  """ transitions are not sorted into bins, returns a vector of all transitions """
  records = get_records('fasta/out.fasta')
  ts = []
  for i in range(0, len(records)):
    ts.extend(count_transitions(labels, records[i]))
  return ts 

