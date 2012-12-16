from collections import Counter

if __name__ == '__main__':
    train = open('train.csv').readlines()
    cnt = Counter()
    cnt.update(t.split(',')[0] for t in train)
    print cnt.most_common(4)
    
    test = open('test.csv').readlines()
    result = ""
    for t in range(len(test)-1):
        result += "2\n"#test[t].split(',')
    print test[0].split(',')
    open('result.csv', 'w').write(result)
    
    result = open('result.csv').readlines()
    