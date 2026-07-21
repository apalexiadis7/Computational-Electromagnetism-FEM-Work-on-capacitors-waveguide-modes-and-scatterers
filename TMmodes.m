function[]=TMmodes(refinements)
% ορισμός μιας απλής κυκλικής γεωμετρίας
gd=[1, 0, 0, 1e-2]';
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
% αρχικοποιούμε με μηδενικά το διάνυσμα των Εz σε όλον τον χώρο. Έχουμε
% οριακές συνθήκες Dirichlet οπότε οι κόμβοι στα όρια είναι άγνωστοι,
% επομένως χρειαζόμαστε τα node_id και index (όρα κώδικα ερωτήματος 2)
node_id=ones(Nn,1);
E0 = zeros(Nn,6);
for i=1:Nd
if e(6,i)==0 || e(7,i)==0
node_id(e(1,i))=0;
node_id(e(2,i))=0;
for j=1:6
E0(e(1,i),j)=0;
E0(e(2,i),j)=0;
end
end
end
% επαναρίθμηση αγνώστων κόμβων με το index
counter=0;
for id=1:Nn
if node_id(id)==1
counter=counter+1;
index(id)=counter;
end
end
figure;
pdeplot(p,e,t);
% ο αριθμός των αγνώστων
Nf=counter;
%αρχικοποίηση των ολικών αραιών πινάκων ακαμψίας και μάζης
S = spalloc(Nf,Nf,7*Nf);
T = spalloc(Nf,Nf,7*Nf);
% Ακολουθεί ο πυρήνας του κώδικα, όπως αναγράφεται στις διαφάνειες
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
Se=zeros(3,3);
%τοπικός πίνακας ακαμψίας
Te=zeros(3,3);
%τοπικός πίνακας μάζας
for i=1:3
for j=1:3
Se(i,j) = (b(i)*b(j)+c(i)*c(j))*Ae;
if i==j
Te(i,j) = Ae/6;
else
Te(i,j) = Ae/12;
end
if (node_id(n(i))==1)
if (node_id(n(j))==1)
S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
T(index(n(i)),index(n(j))) = T(index(n(i)),index(n(j))) + Te(i,j);
end
end
end
end
end
% επίλυση του προβλήματος ιδιοτιμών, λαμβάνοντας τις 10 μικρότερες
sigma='smallestabs';
% με αυτή την έκφραση της συνάρτησης, μας επιστρέφεται ένα διάνυσμα που
% περιέχει τις ιδιοτιμές τις οποίες ακριβώς μετά πλοτάρουμε
d=eigs(S,T,10,sigma);
figure;
plot(d,'o');
% αυτή η έκφραση της συνάρτησης επιστρέφει έναν τετραγωνικό πίνακα D με τις
% ιδιοτιμές του στην κύρια διαγώνιο και έναν πίνακα V με στήλες τα
% αντίστοιχα ιδιοδιανύσματα
[V,D] = eigs(S,T,10,sigma);
% εδώ γίνεται επιλογή των καταλλήλων ιδιοτιμών/ιδιοδιανυσμάτων
V1=zeros(Nf,6);
for i=1:Nf
V1(i,1)=V(i,1);
V1(i,2)=V(i,2);
V1(i,3)=V(i,4);
V1(i,4)=V(i,6);
V1(i,5)=V(i,7);
V1(i,6)=V(i,9);
end
% εδώ "γεμίζουμε" το διάνυσμα Ε0
for i=1:Nn
for j=1:6
if (E0(i,j)==0)&&(node_id(i)==1)
E0(i,j)=V1(index(i),j);
end
end
end
% εδώ παράγονται έξι διαφορετικά διανύσματα, ουσιαστικά οι στήλες του Ε0,
% όπου το καθένα αντιπροσωπέυει έναν συγκεκριμένο ρυθμό
X1=zeros(Nn,1);
X2=zeros(Nn,1);
X3=zeros(Nn,1);
X4=zeros(Nn,1);
X5=zeros(Nn,1);
X6=zeros(Nn,1);
for i=1:Nn
X1(i)=E0(i,1);
X2(i)=E0(i,2);
X3(i)=E0(i,3);
X4(i)=E0(i,4);
X5(i)=E0(i,5);
X6(i)=E0(i,6);
end
% απεικόνιση της z συνιστώσας
for i=1:6
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',E0(1:Nn,i),'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
end
% βρίσκουμε το gradient, εφόσον το εγκάρσιο πεδίο δίνεται από αυτό. Συγκεκριμένα η
% μερική παράγωγος κατά y είναι ανάλογη του Ηx ενώ η μερική παράγωγος κατά
% x έιναι ανάλογη του Ηy
[gradx1,grady1] = pdegrad(p,t,X1);
[gradx2,grady2] = pdegrad(p,t,X2);
[gradx3,grady3] = pdegrad(p,t,X3);
[gradx4,grady4] = pdegrad(p,t,X4);
[gradx5,grady5] = pdegrad(p,t,X5);
[gradx6,grady6] = pdegrad(p,t,X6);
% μας αφορά το μέτρο του πεδίου, το οποίο προκύπτει από τη ρίζα των
% τετραγώνων των Ηx & Ηy
u1=zeros(Ne,1);
u2=zeros(Ne,1);
u3=zeros(Ne,1);
u4=zeros(Ne,1);
u5=zeros(Ne,1);
u6=zeros(Ne,1);
for i=1:Ne
u1(i)=sqrt(((gradx1(i))^2)+((grady1(i))^2));
u2(i)=sqrt(((gradx2(i))^2)+((grady2(i))^2));
u3(i)=sqrt(((gradx3(i))^2)+((grady3(i))^2));
u4(i)=sqrt(((gradx4(i))^2)+((grady4(i))^2));
u5(i)=sqrt(((gradx5(i))^2)+((grady5(i))^2));
u6(i)=sqrt(((gradx6(i))^2)+((grady6(i))^2));
end
% αναπαράσταση των επιθυμητών ρυθμών
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u1,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u2,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u3,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u4,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u5,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
figure;
axis equal;
axis tight;
pdeplot(p,e,t,'XYData',u6,'ColorMap','Jet');
hold on
% να φαίνεται και η κατασκευή
pdegplot(d1);
hold off
end
