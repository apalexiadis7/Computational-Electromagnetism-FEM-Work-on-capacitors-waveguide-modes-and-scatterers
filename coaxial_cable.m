function[capac]=omoaxonikon(refinements)
% μεταβλητή εισόδου ο αριθμός των πυκνώσεων
% μεταβλητή εξόδου η χωρητικότητα
R1=[1, 0, 0, 1.75e-3]';
R2=[1, 0, 0, 0.76e-3]';
gd=[R1,R2];
ns=(char('R1','R2'))';
d1=decsg(gd,'R1-R2',ns);
[p,e,t] = initmesh(d1);
% σε περίπτωση που θέλουμε επιπλέον πύκνωση του πλέγματος
if refinements~=0
for i=1:refinements
[p,e,t]=refinemesh(d1,p,e,t);
end
end
pdeplot(p,e,t);
Nn=size(p,2);
Nd=size(e,2);
Ne=size(t,2);
node_id = ones(Nn,1);
X0 = zeros(Nn,1);
% όλοι οι κόμβοι σε ακμές στον πίνακα e θεωρούνται γνωστοί. Αυτοί που
% βρίσκονται στον εσωτερικό αγωγό είναι υπό τάση 1 Volt ενώ αυτοί του
% εξωτερικού είναι γειωμένοι
for i=1:Nd
node_id(e(1,i))=0;
node_id(e(2,i))=0;
end
% εδώ ο έλεγχος για τις περιοχές με βάση τις συντεταγμένες. Θεωρούμε τον
% μοναδιαίο κύκλο (ακτίνα 1 χιλιοστό) όπου για ακμές εσωτερικές αυτού τις
% θεωρούμε σε δυναμικό 1 volt και εκτός σε 0 volt. Εφόσον ο αρχικοποιημένος
% Χ0 πίνακας έχει μόνον μηδενικά, χρειάζεται να γίνει έλεγχος μόνον για τον
% εσωτερικό αγωγό
for i=1:Nd
if(sqrt((p(1,e(1,i))^2)+(p(2,e(1,i))^2))<1e-3)
X0(e(1,i))=1;
X0(e(2,i))=1;
end
end
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
Se(i,j) = e0*(b(i)*b(j)+c(i)*c(j))*Ae;
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
% επίλυση του συστήματος και μέτρηση του χρόνου
tic
F=S\B;
toc
for i=1:Nn
if ((X0(i)==0)&&(node_id(i)==1))
X0(i)=F(index(i));
end
end
X1=zeros(Nn,1);
%Θέλουμε το gradient του -φ
for i=1:Nn
X1(i)=-X0(i);
end
[gradx,grady] = pdegrad(p,t,X1);
figure;
pdeplot(p,e,t,'XYData',X0,'FlowData',[gradx;grady]);
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
Wn=zeros(Ne,1);
we=0;
% για τον υπολογισμό της χωρητικότητας
for ie = 1:Ne
n(1:3) = t(1:3,ie);
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
Ste=zeros(3,3);
Xt=[X0(n(1)) X0(n(2)) X0(n(3))];
for i=1:3
for j=1:3
Ste(i,j) = e0*(b(i)*b(j)+c(i)*c(j))*Ae;
Wn(ie) = Wn(ie) + X0(n(i))*Ste(i,j)*X0(n(j));
end
end
we=we+Wn(ie);
end
capac=we;
end
