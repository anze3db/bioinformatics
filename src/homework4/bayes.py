import Orange
import orange
import textwrap
import orngEnsemble

def mediana(l):
    l.sort()
    return l[int(len(l)*(1/2.))]
    

if __name__ == '__main__':
    
    
    data = Orange.data.Table("train.csv")
    new_domain = orange.Domain([a for a in data.domain.variables if a.name != 'polII_presence'], data.domain['polII_presence'])
    data = Orange.data.Table(new_domain, data)

    learner = Orange.classification.bayes.NaiveLearner() # orngEnsemble.RandomForestLearner(trees=100, name="forest") 
    classifier = learner(data)
    
    test = Orange.data.Table("test.csv")
    new_domain = orange.Domain([a for a in data.domain.variables], data.domain['polII_presence'])
    test = Orange.data.Table(new_domain, test)
    
    probs = [classifier(inst, orange.GetProbabilities) for inst in test]
    avgs = [sum(p[i] for p in probs)/len(probs) for i in range(3)]
    med  = []
    med.append(mediana([p[0] for p in probs]))
    med.append(mediana([p[1] for p in probs]))
    med.append(mediana([p[2] for p in probs]))
    predictions = []
    for i in range(len(probs)):
        if probs[i][1] >= probs[i][2]:
            predictions.append('2')
        else:
            predictions.append('3')
            
    open('simple_bayes.csv', 'w').write("\n".join(predictions))
    