
a = open("C:/Users/thliao/Desktop/assembly_result.txt").read().split('\n')
a = [_.split('\t')[0] for _ in a if _]
a = a[1:]
with open('C:/Users/thliao/Desktop/list_ids.txt','w') as f1:
    f1.write('\n'.join(a))