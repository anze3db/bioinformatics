def numbers():
    
    to = 10 ** 2

    [i for i in range(1, to)]
    [i for i in range(1, to) if i % 2 == 0]
    [i for i in range(1, to) if i % 3 == 0]
    [i for i in range(1, to) if (i % 3 == 0 and i % 7 == 0)]
    [i for i in range(1, to) if (i % 3 == 0 or  i % 7 == 0)]
    [i ** 2 for i in range(1, to)]
    [i ** 3 for i in range(1, to)]

    print sum(range(1, to))

def sorting():
    
    l = [ ("David", "Naik"), ("Avid", "Naik"), ("Ales", "Stajdohar"), ("Marko", "Erjavec") ]
    l.sort(key=lambda x: (x[1], len(x[0])))
    
    print l


def personalID():
    names = ['David', 'Ales', 'Marko']
    surnames = ['Naik', 'Stajdohar', 'Erjavec']
    births = [1985, 1973, 1999]

    for n,s,b in zip(names,surnames,births):
        print n,s,b


def countwords():
    from collections import Counter
    s = "mama ata teta stric pa lol pa pa pa lol pa kek pa spet mama"
    c = Counter(s.split(' '))
    print c.most_common(2)
    


def local_maximas(l):
    for p,i,n in zip(l[:-2], l[1:-1], l[2:]):
        if p < i and i > n:
            yield i


if __name__ == '__main__':
    numbers()
    sorting()
    personalID()
    countwords()
    print [i for i in local_maximas([1,2,2,3,4,5,4,3,4,3])]