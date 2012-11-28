'''
Created on Jul 15, 2011

@author: smotko
'''

def lexi_permute(str):
    cnt = 0
    while True:
        cnt += 1
        print " ".join(str)
        
        # find j
        n = len(str) - 1 
        j = n - 1
        while str[j] >= str[j+1]:
            j -= 1
            if j == -1:
                return cnt
        
        # increase str j
        l = n
        while str[j] >= str[l]:
            l -= 1
        str = swap(str, j, l)
    
        # reverse
        k = j + 1
        l = n
        while k < l:
            if k < l:
                str = swap(str, k, l)
            k += 1
            l -= 1
        
    
    
def swap(str, i, j):
    if i > j:
        i,j = j,i
    if i == j:
        return str
    return str[:i] + str[j] + str[i+1:j] + str[i] + str[j+1:]    

if __name__ == '__main__':
    
    print lexi_permute("12345")