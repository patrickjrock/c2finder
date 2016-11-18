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
  e = []
  for i in range(0, len(labels)):
    if data[i] != '-':
      state = int(labels[i])
      if state > len(e):
        e.append([])
      e[state-1].append(data[i])
  return e

def emissions(fname):
  records = get_records('fasta/out.fasta')
  es = []
  for i in range(1,len(records)):
    es.append(count_emissions(records[0], records[i]))
  final = []
  for e in es:
    for i, state in enumerate(e):
      if i+1 > len(final):
        final.append([])
      final[i].extend(state)
  return final


def count_transitions(labels, data):
  pass

records = get_records('fasta/out.fasta')
print records
print emissions('fasta.out.fasta')
