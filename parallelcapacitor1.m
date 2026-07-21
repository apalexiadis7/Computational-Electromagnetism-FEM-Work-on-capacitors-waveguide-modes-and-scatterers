function[capac]=
parallelcapacitor1
(wide,thin,dist,volt,eer,refinements,
domainexpansion
)
% με την αναγραφόμενη σειρά οι μεταβλητές εισόδου της συνάρτησης: το φάρδος
% των παραλλήλων πλακών, το πάχος τους, η αναμεταξύ τους απόσταση, η τάση
% (θα υπολογισθεί από +-V/2 για την άνω και κάτω πλάκα αντίστοιχα), η
% διηλεκτρική σταθερά του μέσου πλήρωσης ανάμεσά τους και ο αριθμός των
% πυκνώσεων του πλέγματος (πόσες φορές θα καλέσουμε τη refinemesh). H
% μεταβλητή εξόδου είναι η χωρητικότητα.
% Στις μεταβλητές εισόδου προστέθηκε και μια νέα μεταβλητή η οποία
% καθορίζει το μέγεθος του χωρίου υπολογισμού και χρησιμοποιείται στον
% πίνακα γεωμετρίας
gd=[3 3 3 3;
4 4 4 4;
(-domainexpansion/2)*wide -wide/2 -wide/2 -wide/2;
(domainexpansion/2)*wide wide/2 wide/2 wide/2;
(domainexpansion/2)*wide wide/2 wide/2 wide/2;
(-domainexpansion/2)*wide -wide/2 -wide/2 -wide/2;
(domainexpansion/2)*wide thin+(dist/2) -dist/2 dist/2;
(domainexpansion/2)*wide thin+(dist/2) -dist/2 dist/2;
(-domainexpansion/2)*wide dist/2 -thin-(dist/2) -dist/2;
(-domainexpansion/2)*wide dist/2 -thin-(dist/2) -dist/2];
% ορισμός της γεωμετρίας, οι κόμβοι των τριών περιοχών (άνω & κάτω πλάκα και το ολικό χωρίο)
d1=decsg(gd);
[p,e,t] = initmesh(d1);
% σε περίπτωση που θέλουμε επιπλέον πύκνωση του πλέγματος
if refinements~=0
for i=1:refinements
[p,e,t]=refinemesh(d1,p,e,t);
end
end
Nn=size(p,2);
Nd=size(e,2);
Ne=size(t,2);
node_id = ones(Nn,1);
X0 = zeros(Nn,1);
% έλεγχος περιοχών για άνω και κάτω πλάκα αντίστοιχα
for ie=1:Ne
if(t(4,ie)==2)
node_id(t(1,ie),1)=0;
node_id(t(2,ie),1)=0;
node_id(t(3,ie),1)=0;
X0(t(1,ie))=volt/2;
X0(t(2,ie))=volt/2;
X0(t(3,ie))=volt/2;
elseif(t(4,ie)==3)
node_id(t(1,ie),1)=0;
node_id(t(2,ie),1)=0;
node_id(t(3,ie),1)=0;
X0(t(1,ie))=-volt/2;
X0(t(2,ie))=-volt/2;
X0(t(3,ie))=-volt/2;
end
end
% οπτικοποίηση της κατασκευής χωρίς το πλέγμα αλλά με τις ονομασίες των περιοχών κλπ
figure;
pdegplot(d1,EdgeLabels="on",FaceLabels="on");
% επαναρίθμηση αγνώστων κόμβων με το index
counter=0;
for id=1:Nn
if node_id(id)==1
counter=counter+1;
index(id)=counter;
end
end
% εύρεση των γνωστών κόμβων για μετέπειτα απεικόνισή τους στο πλέγμα
u=zeros(Nn-counter,2);
count=0;
for i=1:Nn
if(node_id(i)==0)
count=count+1;
u(count,1)=p(1,i);
u(count,2)=p(2,i);
end
end
syntx=zeros(length(u));
synty=zeros(length(u));
for i=1:length(u)
syntx(i)=u(i,1);
synty(i)=u(i,2);
end
% Νf ο αριθμός των αγνώστων, ισούται με τον μετρητή του προηγούμενου βρόχου
% ο οποίος αυξανόταν όταν τα node_id ήταν ίσα με τη μονάδα, οπότε επρόκειτο
% για άγνωστο
Nf=counter;
S = spalloc(Nf,Nf,7*Nf);
B = zeros(Nf,1);
e0=8.854e-12;
% διηλεκτρική σταθερά του κενού.
% Ακολουθεί ο πυρήνας του κώδικα, όπως αναγράφεται στις διαφάνειες
for ie = 1:Ne
n(1:3) = t(1:3,ie);
rg = t(4,ie);
x(1:3) = p(1,n(1:3));
y(1:3) = p(2,n(1:3));
De = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
Ae = abs(De/2);
b(1) = (y(2)-y(3))/De;
c(1) = (x(3)-x(2))/De;
b(2) = (y(3)-y(1))/De;
c(2) = (x(1)-x(3))/De;
b(3) = (y(1)-y(2))/De;
c(3) = (x(2)-x(1))/De;
Se=zeros(3,3);
for i=1:3
for j=1:3
% εντός των δύο πλακών έχουμε άλλο διηλεκτρικό, κάνουμε
% έναν επιπλέον έλεγχο όσον αφορά την περιοχή των
% κόμβων
if(rg==4)
Se(i,j) = e0*eer*(b(i)*b(j)+c(i)*c(j))*Ae;
else
Se(i,j) = e0*(b(i)*b(j)+c(i)*c(j))*Ae;
end
if (node_id(n(i))==1)
if (node_id(n(j))==1)
S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
else
B(index(n(i))) = B(index(n(i))) - Se(i,j)*X0(n(j));
end
end
end
end
end
F=mtimes(inv(S),B);
% επίλυση του συστήματος
for i=1:Nn
if X0(i)==0
% για τα άγνωστα στοιχεία, πρέπει να προσθέσουμε τις τιμές τους στο
% διάνυσμα Χ0. Επομένως για κάθε άγνωστο στο αρχικό διάνυσμα Χ0,
% όπου δηλαδή έχουμε μηδενική τιμή, ψάχνουμε το αντίστοιχο στοιχείο
% στο index, όπου υπάρχει η αρίθμηση των αγνώστων βάσει της οποίας
% είναι διατεταγμένα τα μητρώα S & F.
X0(i)=F(index(i));
end
end
% με αυτή την εντολή κατανοούμε και την πύκνωση του πλέγματος οπτικά
figure;
pdeplot(p,e,t);
hold on
% εδώ οπτικοποιούνται και οι γνωστοί κόμβοι που βρήκαμε προηγουμένως
plot(syntx,synty,'LineStyle','none','Marker','o');
hold off
X1=zeros(Nn,1);
%Θέλουμε το gradient του -φ
for i=1:Nn
X1(i)=-X0(i);
end
[gradx,grady] = pdegrad(p,t,X1);
figure;
pdeplot(p,e,t,'XYData',X0,'FlowData',[gradx;grady],'ColorMap','Parula');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
%ένα πιο πυκνό διάγραμμα του πεδίου μόνον
Fu=pdeInterpolant(p,t,X1);
[X,Y]=meshgrid((-5/2)*wide:(5/100)*wide:(5/2)*wide);
u=evaluate(Fu,X,Y);
u_shaped=reshape(u,size(X));
[gx,gy]=gradient(u_shaped);
figure;
quiver(X,Y,gx,gy);
hold on
pdegplot(d1);
hold off
% όσον αφορά και τον υπολογισμό της χωρητικότητας τώρα, θα εκμεταλλευτούμε
% τον τύπο W=(1/2)*C*V^2
work=zeros(1,Ne);
% μητρώο της ενέργειας
for i=1:Ne
% Τα στοιχεία του μητρώου θα
% είναι τα τετράγωνα των τιμών του ηλεκτρικού πεδίου στην επιθυμητή
% περιοχή, τα οποία θα πάρουμε από το [gradx,grady]
n(1:3) = t(1:3,i);
rg = t(4,i);
x(1:3) = p(1,n(1:3));
y(1:3) = p(2,n(1:3));
De = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
Ae = abs(De/2);
if rg==4
work(1,i)=e0*eer*(((gradx(i))^2)+((grady(i))^2))*Ae;
else
work(1,i)=0;
end
end
% αθροίζουμε τις τιμές του πίνακα work για την διακριτή ολοκλήρωση
we=0;
for i=1:Ne
we=we+work(1,i);
end
capac=we/((volt)^2);
end
