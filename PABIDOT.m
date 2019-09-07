clear all

sdev=0.3;%Standard deviation of the noise added.

[dataset]=readtable('input_data.csv'); %Path to the input dataset

data=dataset{:,1:size(dataset,2)-1};
classrange=dataset{:,size(dataset,2)};
path = ['PABIDOT_perturbed.csv']; % Path to the perturbed (output) dataset


gpastdev=std(data);
gpamean=mean(data);

c=zscore(data);

comatoriginal=cov(c);

attcount=size(c,2);
comsincos=nchoosek(1:attcount,2); %Total sin cos replacements

optimalprivacy=0;
for tempaxis=1:attcount
    randaxis=tempaxis;
    identityref=eye(attcount);
    identityref(:,randaxis)=-identityref(:,randaxis);
    
    %%%%%%%%%%%%%%Evaluating theta for 1:180 degrees%%%%%%%%%%%%%%%%%%%%%%%%%%
    for rotangle=1:180
        if(rotangle==30||rotangle==45||rotangle==60||rotangle==90||rotangle==120||rotangle==135||rotangle==150||rotangle==180)
           continue 
        end
        %%%%%%%%%%%%%%%%%%%%%%Generating the rotational matrix%%%%%%%%%%%%%%%%%%%%%    
        radianvalue=rotangle*pi/180;
        anglerad=radianvalue;
      
        identitymatinit=abs(eye(attcount));
        for i=1:size(comsincos,1)%Brows through all the possible combinations (matrices)
           coordinates=combvec(comsincos(i,:),comsincos(i,:))';%coordinates for sin cos
           idenitymatTemp=abs(eye(attcount));
           idenitymatTemp(coordinates(1,1),coordinates(1,2))=cos(anglerad);
           idenitymatTemp(coordinates(2,1),coordinates(2,2))=sin(anglerad);
           idenitymatTemp(coordinates(3,1),coordinates(3,2))=-sin(anglerad);
           idenitymatTemp(coordinates(4,1),coordinates(4,2))=cos(anglerad);
           identitymatinit=identitymatinit*idenitymatTemp;   
        end
        sepval(rotangle,tempaxis)=min(1+(diag(identitymatinit*(identityref*comatoriginal*identityref')*identitymatinit'))-2*(sum(((comatoriginal*identityref).*identitymatinit)')'));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%%%%%%%%%%%%%%Retrieving optimal theta and optimal axis of reflection%%%%%%%%
[M,I] = min(sepval');
[maxval besttheta]=max(M);
refaxis=I(besttheta);

%%%%%%%%%%%%%%%%Generating the reflection matrix%%%%%%%%%%%%%%%%%%%%%%%%%%
randaxis=refaxis;
identityref=eye(attcount);
identityref(:,randaxis)=-identityref(:,randaxis);

%%%%%%%%%%%%%%Applying reflection transformation to data%%%%%%%%%%%%%%%%%%
bo=identityref*c';
bo=bo';

%%%%%%%%%%%%%%%%Translation Matrix Generation/Application%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
identitytrans=eye(attcount+1); %attcount to add a new row with the homogeneous coordinate
randomnoise=1*rand((attcount),1);%Uniform random nose
identitytrans(1:attcount,(attcount+1))=randomnoise;%Uniform Noise last column
multitrans=bo;
multitrans=horzcat(multitrans,ones(1,size(bo,1))'); %adding an additional column of zeros, inorder to add a new line for homogeneous coordinate
transresult=identitytrans*multitrans';
transresult=transresult';
multitrans=transresult(1:size(bo,1),1:size(bo,2));%removing homogeneous coordinate
bo=multitrans;

%final=identitymatinit*b';
radianvalue2=besttheta*pi/180;
anglerad2=radianvalue2;
identitymatinit2=abs(eye(attcount));
for i=1:size(comsincos,1)%Brows through all the possible combinations (matrices)
   coordinates=combvec(comsincos(i,:),comsincos(i,:))';%coordinates for sin cos
   idenitymatTemp2=abs(eye(attcount));
   idenitymatTemp2(coordinates(1,1),coordinates(1,2))=cos(anglerad2);
   idenitymatTemp2(coordinates(2,1),coordinates(2,2))=sin(anglerad2);
   idenitymatTemp2(coordinates(3,1),coordinates(3,2))=-sin(anglerad2);
   idenitymatTemp2(coordinates(4,1),coordinates(4,2))=cos(anglerad2);
   identitymatinit2=identitymatinit2*idenitymatTemp2;   
end

%%%%%%%%%%%%%%%%%%%%%%%%Application of rotational transformation%%%%%%%
bo=identitymatinit2*bo';

%randomized expansion
s=sign(bo);
bo=abs(bo)+abs((randn(size(bo,1),size(bo,2)).*sdev));
bo=bo.*s;
bo=bo';

bo=bo.*gpastdev+gpamean; %inverse normlization on the dataset

%%%%%%%%%%%%%%%%%%%%%%Applying Stress over the perturbed data %%%%%%%%%
permset=randperm(size(bo,1));
shuffledset=bo(permset,:);
shuffledclasses=classrange(permset,:);

for i = 1:size(bo,2)
   c=['a',num2str(i)];
   att(i)=string(c);
end


fid = fopen(path, 'w') ;
fprintf(fid, '%s,', att) ;
fprintf(fid, 'class\n', att); 
for j=1:size(bo,1)
    fprintf(fid, '%f,', shuffledset(j,:)) ;
    fprintf(fid, '%s\n', cell2mat(shuffledclasses(j))) ;
end
fclose(fid);

