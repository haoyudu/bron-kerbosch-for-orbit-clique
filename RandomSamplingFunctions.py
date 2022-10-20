import random

def RandomSkewSet(Option,SkewSet,Print="False"):
    OriginalLength = len(SkewSet)
    while not(Option=="minimized") and not(Option=="maximized") and not(Option=="random"):
        print("Please input one of the three valid options: minimized, maximized, or random.")
        Option = input();
    RemainingPossibleLines = CommonNeighbours(SkewSet,graph);
    if len(RemainingPossibleLines)==0:
        return SkewSet
    if Option == "random":
        start_time = time.time()
        while len(RemainingPossibleLines)>0:
            counter = len(RemainingPossibleLines)
            Index = random.randint(0,len(RemainingPossibleLines)-1);
            SkewSet.append(RemainingPossibleLines[Index]);
            RemainingPossibleLines.remove(RemainingPossibleLines[Index])
            RemainingPossibleLines = [x for x in RemainingPossibleLines if x in Neighbours[SkewSet[-1]]];
    elif Option == "maximized" or Option == "minimized":
        start_time = time.time()
        Index =random.randint(0,len(Lines)-1);
        '''
        for i in range(2):
            print(len(SkewSet),len(RemainingPossibleLines),round(time.time() - start_time,2));
            start_time = time.time();
            Index =random.randint(0,len(RemainingPossibleLines)-1);
            SkewSet.append(RemainingPossibleLines[Index]);
            RemainingPossibleLines.remove(RemainingPossibleLines[Index])
            NewRemainingPossibleLines = [line for line in RemainingPossibleLines if line in Neighbours(Index,graph)];
            RemainingPossibleLines = NewRemainingPossibleLines;
        '''
        while len(RemainingPossibleLines)>0:
            counter = len(RemainingPossibleLines)
            #print(len(SkewSet),len(RemainingPossibleLines),round(time.time() - start_time,2));
            start_time = time.time();
            Chance = random.randint(0,5);
            if Chance<50:
                Min = len(RemainingPossibleLines);
                Max = 0;
                PossibleNextSkewLines = [];
                for i in range(0,len(RemainingPossibleLines)):
                    Temp = len([line for line in RemainingPossibleLines if line in Neighbours[RemainingPossibleLines[i]]]);
                    if (Temp < Min and Option =="minimized") or (Temp > Max and Option =="maximized"):
                        Min = Temp;
                        PossibleNextSkewLines = [RemainingPossibleLines[i]];
                    elif (Temp == Min and Option =="minimized") or (Temp == Max and Option =="maximized"):
                        PossibleNextSkewLines.append(RemainingPossibleLines[i]);
                Index = random.randint(0,len(PossibleNextSkewLines)-1);
                SkewSet.append(PossibleNextSkewLines[Index]);
                RemainingPossibleLines.remove(PossibleNextSkewLines[Index])
                RemainingPossibleLines = [line for line in RemainingPossibleLines if line in Neighbours[PossibleNextSkewLines[Index]]];
            else:
                counter = len(RemainingPossibleLines)
                Index = random.randint(0,len(RemainingPossibleLines)-1);
                SkewSet.append(RemainingPossibleLines[Index]);
                RemainingPossibleLines.remove(RemainingPossibleLines[Index])
                RemainingPossibleLines = [x for x in RemainingPossibleLines if x in Neighbours[SkewSet[-1]]];        
    if not(Print):
        if q ==4:
            if (len(SkewSet)<13  or len(SkewSet)>25) and not(OriginalLength ==len(SkewSet)):
                print(len(SkewSet),sorted(SkewSet))
        elif q ==5:
            if (len(SkewSet)<17  or len(SkewSet)>47) and not(OriginalLength ==len(SkewSet)):
                print(len(SkewSet),sorted(SkewSet))
        elif q ==7:
            if (len(SkewSet)<37  or len(SkewSet)>88) and not(OriginalLength ==len(SkewSet)):
                print(len(SkewSet),sorted(SkewSet))
        elif q ==8:
            if (len(SkewSet)<46  or len(SkewSet)>99) and not(OriginalLength ==len(SkewSet)):
                print(len(SkewSet),sorted(SkewSet))
        elif Print == True or len(SkewSet) >(3/2)*q**2-(1/2)*q+1 or len(SkewSet)<q**2+1 or (p==2 and len(SkewSet) >q**2+q):
            print(len(SkewSet),sorted(SkewSet))
    return SkewSet;

# This finds a fixed number of random skew sets and outputs how many of each size occur.
# number is a positive integer
def RandomSkewSetLengths(number,Option="random"):
    start_time = time.time()
    Lengths = [];
    NumberOfLengths = 0;
    for i in range(0,number):
        Lengths.append(len(RandomSkewSet(Option,[0,q+2,2*q+4])));
        if i%100000 == 99999:
            print(str(i+1)+" of "+str(number)+" has taken "+str(round(time.time() - start_time,2))+" seconds.")
    for k in range(1,max(Lengths)+1):
        x = 0;
        for i in range(0,number):
            if Lengths[i]==k:
                x = x+1;
        if not(x==0):
            NumberOfLengths = NumberOfLengths+x;
            print("There are "+str(x)+" skew sets of size "+str(k)+" in this collection.")
    if not(NumberOfLengths == number):
        print("An error occured");

def MinimalSkewSetMCMC(NumberOfIterations,SkewSet=RandomSkewSet("random",[0,q+2,2*q+4])):
    SkewSetLengths = [len(SkewSet)]
    iteration = 0;
    #counter = 0;
    #DeadLoops = [0]
    while iteration < NumberOfIterations:
        ShortNewSkewSet = [];
        for index in range(3,len(SkewSet)):
            Chance = random.randint(0,3);
            if Chance ==0:
                ShortNewSkewSet.append(SkewSet[index])
        NewSkewSet = RandomSkewSet("random",ShortNewSkewSet)
        if not(len(NewSkewSet)>len(SkewSet)):
            if len(NewSkewSet)<len(SkewSet):
                print("Size of Skew Set:", len(NewSkewSet),"Iteration:",iteration)
            SkewSet=NewSkewSet
            SkewSetLengths.append(len(SkewSet))
            iteration = iteration+1
        #    DeadLoops.append(0)
        #else:
        #    DeadLoops[-1]=DeadLoops[-1]+1
        if (q == 2 and len(SkewSet) == 5) or (q == 3 and len(SkewSet) ==7) or (q == 4 and len(SkewSet) == 13):
            print("Iterations: ", iteration)
            iteration = NumberOfIterations
    print(counter)    
    if len(SkewSetLengths)>100:
        print([SkewSetLengths[y] for y in range(-100,0)])
    else:
        print(SkewSetLengths)
    #if len(DeadLoops)>100:
    #    print([DeadLoops[y] for y in range(-100,0)])
    #else:
    #    print(DeadLoops)
    for x in range(min(SkewSetLengths),max(SkewSetLengths)+1):
        counter = 0;
        for i in range(len(SkewSetLengths)):
            if SkewSetLengths[i]==x:
                counter = counter+1;
        print([x,counter])
    print(sorted(SkewSet))


def MinimalSkewSetMCMC(k,Option,NumberOfIterations,SkewSet=RandomSkewSet("random",[0,q+2,2*q+4])):
    SkewSetLengths = [len(SkewSet)]
    iteration = 0;
    while iteration < NumberOfIterations:
        if iteration % 1000 == 999:
            iteration = iteration+1
            print("Iteration:", iteration)
        ShortNewSkewSet = [];
        for index in range(3,len(SkewSet)):
            Chance = random.randint(0,k);
            if not(Chance ==0):
                ShortNewSkewSet.append(SkewSet[index])
        NewSkewSet = RandomSkewSet(Option,ShortNewSkewSet)
        if not(sorted(NewSkewSet) == sorted(SkewSet)) and len(NewSkewSet)<20 and q==5:
            print(len(NewSkewSet),sorted(NewSkewSet))
        if not(len(NewSkewSet)>len(SkewSet)) and not(sorted(NewSkewSet) == sorted(SkewSet)):
            if len(NewSkewSet)<len(SkewSet):
                print("Size of Skew Set:", len(NewSkewSet),"Iteration:",iteration)
            SkewSet=NewSkewSet
            SkewSetLengths.append(len(SkewSet))
            iteration = iteration+1
        if (q == 2 and len(SkewSet) == 5) or (q == 3 and len(SkewSet) ==7):
            print("Iterations: ", iteration)
            iteration = NumberOfIterations
    for x in range(min(SkewSetLengths),max(SkewSetLengths)+1):
        counter = 0;
        for i in range(len(SkewSetLengths)):
            if SkewSetLengths[i]==x:
                counter = counter+1;
        print([x,counter])
    print(sorted(SkewSet))


def MaximalSkewSetMCMC(k,Option,NumberOfIterations,SkewSet=RandomSkewSet("random",[0,q+2,2*q+4])):
    SkewSetLengths = [len(SkewSet)]
    iteration = 0;
    while iteration < NumberOfIterations:
        if iteration % 1000 == 999:
            iteration = iteration+1
            print("Iteration:", iteration)
        ShortNewSkewSet = [];
        for index in range(3,len(SkewSet)):
            Chance = random.randint(0,k);
            if not(Chance ==0):
                ShortNewSkewSet.append(SkewSet[index])
        NewSkewSet = RandomSkewSet(Option,ShortNewSkewSet)
        if not(sorted(NewSkewSet) == sorted(SkewSet)) and len(NewSkewSet)<20 and q==5:
            print(len(NewSkewSet),sorted(NewSkewSet))
        if not(len(NewSkewSet)<len(SkewSet)) and not(sorted(NewSkewSet) == sorted(SkewSet)):
            if len(NewSkewSet)>len(SkewSet):
                print("Size of Skew Set:", len(NewSkewSet),"Iteration:",iteration)
            SkewSet=NewSkewSet
            SkewSetLengths.append(len(SkewSet))
            iteration = iteration+1
        if (q == 2 and len(SkewSet) == 5) or (q == 3 and len(SkewSet) ==7):
            print("Iterations: ", iteration)
            iteration = NumberOfIterations
    for x in range(min(SkewSetLengths),max(SkewSetLengths)+1):
        counter = 0;
        for i in range(len(SkewSetLengths)):
            if SkewSetLengths[i]==x:
                counter = counter+1;
        print([x,counter])
    print(sorted(SkewSet))
