import csv

def read_genes(file_path):
    file=open(file_path,"r")
    gene_dict={}
    header=""
    sequence=""
    for line in file:
        #Checked if the line is a header or sequnce and made decisions according to that.
        if line.startswith(">"):
            gene_dict[header.strip("\n")]=sequence.strip("\n")
            header=line
            sequence=""
        else:
            sequence+=line
    gene_dict[header.strip("\n")]=sequence.strip("\n")
    file.close()
    #Got rid of the item which we don't want it to be in our dictonary.
    gene_dict.pop("")
    return gene_dict

def get_fragments(gene_dict, frag_len = 50):
    frag_dict = {}
    for item in gene_dict:
        #We split and format our header to get integer values.
        temp=item.split("|")
        itemHeaderchar=temp[0]
        temp=str(temp[1].strip())
        temp=temp.split("-")
        diff=int(temp[1])-int(temp[0])
        #Variable 'value' is our sequence.
        value=gene_dict[item]
        #With using this while loop we filtered out the one's which's lenghts are smaller than 50.
        while diff >= frag_len:
            tempStr = ""
            #We get first 50(frag_len) character.
            for i in range(frag_len):
                tempStr += value[i]
            #We delete the first 50(frag_len) character to keep dividing them.
            value=value[frag_len:]
            #We add 50(frag_len) to our header's numbers.
            temp[0]=int(temp[0])+frag_len
            #We reformat the header to save with right values.
            itemHeader=itemHeaderchar+"|"+str(temp[0]-frag_len)+"-"+str(temp[0])
            frag_dict[itemHeader]=tempStr
            #We change the difference to see if we have at least 50 characters and prevent the while to be an infinite loop.
            diff=int(temp[1])-int(temp[0])
    return frag_dict
                
def filter_frags(frag_dict, threshold = 0.7):
    def get_similarity(s1, s2):
        counter = 0
        #Just to be sure :)
        if len(s1) == len(s2):
            for i in range(len(s1)):
                #To check how many characters are same.
                if s1[i]==s2[i]:
                    counter+=1
        return (counter/len(s1))

    value_list = []
    key_list = []
    temp_value_list = []
    dissimilar_frag_dict = {}
    is_similar = 0

    #We first convert the dictonary to 2 lists which are for keys and values.
    for frag_key in frag_dict:
        value_list.append(frag_dict[frag_key])
        key_list.append(frag_key)

    #Checked for the similar ones.
    for i in range(len(value_list)):
        #Our range starts from i to len(value_list) to match the result in template.
        for j in range(i,len(value_list)):
            if get_similarity(str(value_list[i]),str(value_list[j])) >= threshold:
                is_similar += 1
        #To see if they are similar we needed to check for at least 2 values. One for itself, one for the similar one which is different.
        if is_similar < 2:
            temp_value_list.append(value_list[i])
            dissimilar_frag_dict[key_list[i]] = value_list[i]
        is_similar = 0

    return dissimilar_frag_dict

def get_sentences(dissimilar_frag_dict):
    def generate_kmers(seq, k):
        result_string = ""
        while len(seq) >= k:
            #We add first 4(k) characters to a string.
            for i in range(k):
                result_string += seq[i]
            #We delete the first one to keep getting the first 4(k) characters.
            seq=seq[1:]
            result_string+=" "
        return result_string.strip()
    
    sentences_dict = {}
    #For every value we generate k_mers and equalize it to its key in a new dictonary.
    for item in dissimilar_frag_dict:
        sentences_dict[item] = (generate_kmers(dissimilar_frag_dict[item],4))
    return sentences_dict

def clean_dict(sentences_dict):
    def clean_sentence(sentence):
        #We first get every word of a sentence and put them in a list.
        list_of_sentence=sentence.split(' ')
        temp_sentence = []
        result_sentence = ""
        for i in list_of_sentence:
            #If we don't have the word in out list, we add it.
            if i not in temp_sentence:
                temp_sentence.append(i)
        #We add them to a string again to return it as the result.
        for i in temp_sentence:
            result_sentence+=i+" "
        return result_sentence.strip()
    clean_sentences_dict = {}
    #For every we clean the sentence and equalize it to its key in a new dictonary.
    for item in sentences_dict:
        clean_sentences_dict[item]=clean_sentence(sentences_dict[item])
    return clean_sentences_dict

def write_genes(file_path, clean_sentences_dict):
    file = open(file_path,"w",newline="")
    writer=csv.writer(file)
    #This is our header to categorize the cells
    header = ['fragment_id', 'sentence' ,'sentence_length' ,'number_of_words']
    writer.writerow(header)
    #For every key-value pair, we write it to our file.
    for item in clean_sentences_dict:
        data = [item, clean_sentences_dict[item], len(clean_sentences_dict[item]), len(clean_sentences_dict[item].split( ))]
        writer.writerow(data)
    file.close()
    

def main():
    #These are all to write the results that we want to see(from template).
    genes=read_genes("input.txt")
    print(len(genes))
    fragments=get_fragments(genes)
    print(len(fragments))
    filtered=filter_frags(fragments)
    print(len(filtered))
    sentences=get_sentences(filtered)
    word_count = len(sentences[next(iter(sentences))].split(" "))
    print(word_count)
    first_sentences_lenght = 0
    for item in sentences:
        first_sentences_lenght += len(sentences[item].split(" "))
    cleaned=clean_dict(sentences)
    last_sentence_length = 0
    for item in cleaned:
        last_sentence_length += len(cleaned[item].split(" "))
    print(first_sentences_lenght-last_sentence_length)
    write_genes("output.csv",cleaned)
        

#And finally we needed to call the main for our code to work :)a
main()
