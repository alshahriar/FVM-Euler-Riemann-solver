% In this brach it will be developed for HLLC

%% Preprocessing
clear all
% Generating MATLAB Decomposed Geometry
up=5.5; down = - 5.5; left=-6; right=7.5;
%upper,down and exit boundary
%[code numberOfPoints x-coordinates y-coordinates]
R1 = [3;4; left;right;right;left; up;up;down;down];
%circle(large) at the inlet boundary
%[code x-center y-center redius]
C1 = [1; 3.5; 0; 8];
%circle(small) cylinder wall boundary
C2 =[1; 0; 0; 1];
P1=[2;5; -1;5;5;0;-1; 4;4;5/sqrt(3);0;0];
%[equalizing the dimentions]
C1 = [C1;zeros(length(R1) - length(C1),1)];
C2 = [C2;zeros(length(R1) - length(C2),1)];
gd = [R1,C1,C2];
ns = char('R1','C1','C2'); ns = ns';
sf = 'R1-(R1-C1)-C2';
%defining decomposed geometry
[domain] = decsg(gd,sf,ns); pdegplot(domain); axis('equal');
%initiate mesh and mesh refining
[coordinate,initEdge,tri] = initmesh(domain);
[coordinate,initEdge,tri]=refinemesh(domain,coordinate,initEdge,tri);
[coordinate,initEdge,tri]=refinemesh(domain,coordinate,initEdge,tri);
[coordinate,initEdge,tri]=refinemesh(domain,coordinate,initEdge,tri);
figure(1); pdeplot(coordinate,initEdge,tri); axis('equal');
%pause(0.8); close(figure(1));
%domain and meshing information
nV = max(size(coordinate)); % total number of grid points
nBE = max(size(initEdge)); % number of boundary edges
nT = max(size(tri)); % number of Trianguler mesh
tInfo = tri(1:3,:);
nE = (3*nT - nBE)/2; %number of internal edges
nTotalE=nE+nBE;
edge = zeros(8,nE);%edge lists
pData = zeros(nV,nV);
edge(1:2,nE+1:nTotalE)=initEdge(1:2,:);
edge(5,1:nE) = 0;
% finding edges that constitute the cylinder
distance = sqrt(coordinate(1,edge(1,nE+1:nTotalE)).^2 + ...
coordinate(2,edge(1,nE+1:nTotalE)).^2);
edge(5,nE+1:nTotalE) = (distance < 1.1); % 1 for wall

for i = nE+1:nTotalE,
xyT(:,:)=[coordinate(:,edge(1,i)),coordinate(:,edge(2,i))];
    if(xyT(1,1) == xyT(1,2)),
    edge(5,i)=2; % 2 for exit
    line(xyT(1,:),xyT(2,:),'color','y');
    end
    if(edge(5,i)==1), 
        line(xyT(1,:),xyT(2,:),'color','r');
    end
    if(edge(5,i)==0), 
        line(xyT(1,:),xyT(2,:),'color','g');
    end
end

% setting data for interior edges
for indexBoundary = nE+1:nTotalE,
firstPoint = edge(1,indexBoundary);
secondPoint = edge(2,indexBoundary);
% finding min and max index of points
minIndex = min(firstPoint, secondPoint);
maxIndex = max(firstPoint, secondPoint);
% Adding edge to list noting triangle
pData(minIndex,maxIndex) = -indexBoundary;
end
internalEdge = 0;
% browsing over all triangles and find all edges
for indexTri = 1:nT,
    for j = 1:3,
        firstPoint = tInfo(j, indexTri);
        secondPoint = tInfo(mod(j,3)+1, indexTri);
        minIndex = min(firstPoint, secondPoint);
        maxIndex = max(firstPoint, secondPoint);
        tIndex = pData(minIndex,maxIndex);
        if (tIndex == 0)
            pData(minIndex,maxIndex) = indexTri;
        elseif (tIndex > 0),
            internalEdge = internalEdge + 1;
            edge(1,internalEdge) = minIndex;
            edge(2,internalEdge) = maxIndex;
            if (minIndex == firstPoint)
                edge(3,internalEdge) = indexTri;
                edge(4,internalEdge) = tIndex;
            else
            % indexTri is the index of the triangle situated right
            % side of the edge
                edge(3,internalEdge) = tIndex;
                edge(4,internalEdge) = indexTri;
            end
            pData(minIndex,maxIndex) = 0;
        else
            edge(1,-tIndex) = firstPoint;
            edge(2,-tIndex) = secondPoint;
            edge(3,-tIndex) = indexTri;
            pData(minIndex,maxIndex) = 0;
        end
    end
end
clear pData minIndex maxIndex internalEdge tIndex rb firstPoint;
clear secondPoint e t indexTri domain j R1 C1 C2 ns gd sf;
%finding normal derictions
edge(6:7,:)=[coordinate(2,edge(1,:))-coordinate(2,edge(2,:)); ...
coordinate(1,edge(2,:))-coordinate(1,edge(1,:))];
edge(8,:) = sqrt(edge(6,:).^2 + edge(7,:).^2);
edge(6:7,:) = [edge(6,:)./edge(8,:); edge(7,:)./edge(8,:)];
% Flow field properties
MachFree = 1.75; % freestream mach number
MachInit = 0.75; % initial value - Mach
gamma = 1.4;
CD = 0; CL = 0; % inital lift and drag
timeStep=zeros(1,nT); % varible for local time stepping
% freestream state vector
UFreeStream = [1; MachFree; 0; 1/(gamma-1)/gamma + 0.5*MachFree^2];
%initial flow field data
U(1,1:nT) = 1;
U(2,1:nT) = MachInit;
U(3,1:nT) = 0;
U(4,1:nT) = 1/(gamma-1)/gamma + 0.5*MachInit^2;
CFL = 0.9; n=0;
FBalance=zeros(4,nT); %Flux Balance
nIteration = 300;



%% performing iteration
while (n<nIteration)
% browsing over the all the edges
    for i = 1:nTotalE,
        iL = edge(3,i);
        if (i > nE) % boundary edges
            if (edge(5,i) == 0), % farfield
                [F, maxSpeed] = HLin(edge(8,i),edge(6:7,i),U(:,iL),UFreeStream);
            elseif(edge(5,i) == 1) % circle
                pVel=sqrt(U(2,iL)^2 + U(3,iL)^2)/U(1,iL);
                normalVel = edge(6,i)*(U(2,iL)/U(1,iL)) + edge(7,i)*(U(3,iL)/U(1,iL));
                tengentVel = sqrt(pVel^2 - normalVel^2);
                pressure = (gamma-1)*(U(4,iL) - 0.5*U(1,iL)*tengentVel^2);
                aW = sqrt(gamma*pressure/(U(1,iL)));
                maxSpeed = edge(8,i)*(abs(normalVel) + aW);
                F = edge(8,i)*[0; pressure*edge(6,i); pressure*edge(7,i); 0]; 
            elseif (edge(5,i) == 2), % exit
                [F, maxSpeed] = HLin(edge(8,i),edge(6:7,i),U(:,iL),U(:,iL)); 
            end
                FBalance(:,iL) = FBalance(:,iL) - F; timeStep(iL) = timeStep(iL) + maxSpeed;
        else % internal edges 
            iR = edge(4,i);
            [F, maxSpeed] = HLin(edge(8,i),edge(6:7,i),U(:,iL),U(:,iR)); 
            FBalance(:,iL) = FBalance(:,iL) - F;
            FBalance(:,iR) = FBalance(:,iR) + F; 
            timeStep(iL) = timeStep(iL) + maxSpeed; 
            timeStep(iR) = timeStep(iR) + maxSpeed;
        end
    end
%time step
dtA = 2*CFL./timeStep;  %dtA = 0.2;

%Updating the current solution  
U(1,:) = U(1,:) - dtA.*FBalance(1,:);
U(2,:) = U(2,:) - dtA.*FBalance(2,:);
U(3,:) = U(3,:) - dtA.*FBalance(3,:);
U(4,:) = U(4,:) - dtA.*FBalance(4,:);

FBalance=FBalance*0;  timeStep=timeStep*0; n=n+1;
end
%% Post Processing starts

velocityVector = sqrt(U(2,:).^2 + U(3,:).^2)./U(1,:);
Pr = (gamma-1)*(U(4,:) - 0.5*U(1,:).*velocityVector.^2); a = sqrt(gamma*Pr./U(1,:));
Mach = velocityVector./a;
%Drag and lift
CDrag = 0; CLift = 0; for i = nE+1:nTotalE,
iL = edge(3,i);
if (edge(5,i) == 1), % around the circle 
pVel=sqrt(U(2,iL)^2 + U(3,iL)^2)/U(1,iL);
normalVel = edge(6,i)*(U(2,iL)/U(1,iL)) + edge(7,i)*(U(3,iL)/U(1,iL)); tengentVel = sqrt(pVel^2 - normalVel^2);
pressure = (gamma-1)*(U(4,iL) - 0.5*U(1,iL)*tengentVel^2); aW = sqrt(gamma*pressure/(U(1,iL)));
maxSpeed = edge(8,i)*(abs(normalVel) + aW);
F = edge(8,i)*[0; pressure*edge(6,i); pressure*edge(7,i); 0]; CDrag = CDrag - F(2); CLift = CLift - F(3);
end
end
CDrag = CDrag/(MachFree^2); CLift = CLift/(MachFree^2);

%Plotting
figure; clf; p2=patch('Vertices',coordinate','Faces',tInfo','CData',Mach'); 
set(p2,'FaceColor','flat','EdgeColor','none','LineWidth',0.5) 
title(sprintf('Mach Profile for Free Stream Velocity of %5.2f, Drag =%5.2f',...
MachFree,CDrag));
axis([left right down up]); colorbar; %caxis([0.15,0.4]);

figure; clf; p2=patch('Vertices',coordinate','Faces',tInfo','CData',Pr'); 
set(p2,'FaceColor','flat','EdgeColor','none','LineWidth',0.5)
title(sprintf('Pressure Profile for Free Stream Velocity of %5.2f,Drag=%5.2f',...
MachFree,CDrag));
axis([left right down up]); colorbar; %caxis([0.15,0.4]);
% Vector Plot 
centerP=zeros(2,nT);
xx=[0 0 0];    yy=[0 0 0];
for i=1:nT
for j = 1:3
xx(j) = coordinate(1,tInfo(j,i));
yy(j) = coordinate(2,tInfo(j,i));
end
centerP(1,i)= mean(xx);
centerP(2,i)= mean(yy);
end
x2=randi(nT,[1,nT/2]);
figure;
quiver(centerP(1,x2(:)),centerP(2,x2(:)),U(2,x2(:)),U(3,x2(:)),1.5,'k');
axis([left right down up]);
for i = nE+1:nTotalE,
xyT(:,:)=[coordinate(:,edge(1,i)),coordinate(:,edge(2,i))];
line(xyT(1,:),xyT(2,:),'color','b');
end