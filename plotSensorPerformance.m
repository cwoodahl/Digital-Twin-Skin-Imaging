%This code can be used to pull raw pixel data or image data from a sensor
%or processed imaged

sensorChoice = sensorGet(sensor,'name');
dataIPcurr = ipGet(ip,'data display');

[r c] = size(dataIPcurr);

%pull and plot processed image values
epi1DataR = sum(dataIPcurr(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
epi1DataG = sum(dataIPcurr(1+1:round(r/2)-1,:,2))/length(1+1:round(r/2)-1);
epi1DataB = sum(dataIPcurr(1+1:round(r/2)-1,:,3))/length(1+1:round(r/2)-1);

epi6DataR = sum(dataIPcurr(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
epi6DataG = sum(dataIPcurr(round(r/2)+2:r-1,:,2))/length(round(r/2)+2:r-1);
epi6DataB = sum(dataIPcurr(round(r/2)+2:r-1,:,3))/length(round(r/2)+2:r-1);

figure(99);
hold on;
plot(epi1DataG,'LineWidth',2);
plot(epi6DataG,'LineWidth',2);
title('Average Green Pixel Value, Type I Epidermis');
legend({'Large Blood Volume','Small Blood Volume'});

figure(100);
hold on;
%plot(epi1DataR,'LineWidth',2);
plot(epi1DataR,'LineWidth',2,'color','red');
plot(epi1DataG,'LineWidth',2,'color','green');
plot(epi1DataB,'LineWidth',2,'color','blue');
title('Epidermis Type I Sensor Averages')



figure(101);
hold on;
%plot(epi6DataR,'LineWidth',2);
plot(epi6DataR,'LineWidth',2,'color','red');
plot(epi6DataG,'LineWidth',2,'color','green');
plot(epi6DataB,'LineWidth',2,'color','blue');
title('Epidermis Type VI Sensor Averages');

%depending on sensor used put a 1 in if statement to run
%monochrome
if(0)
    senDat = sensor.data.volts;
    [r c] = size(senDat);
    epi1Data = sum(senDat(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
    epi6Data = sum(senDat(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);

    figure; hold on;
    plot(epi1Data,'LineWidth',2);
    plot(epi6Data,'LineWidth',2);
end

%imx363
if(0)
    senDat = sensor.data.volts;
    [r c] = size(senDat);
    senDatR = senDat(1:2:r,1:2:c);
    senDatB = senDat(2:2:r,2:2:c);
    senDatG = senDat(1:2:r,2:2:c) + senDat(2:2:r,1:2:c);
    figure; imagesc(senDatR); colorbar;
    figure; imagesc(senDatB); colorbar;
    figure; imagesc(senDatG); colorbar;
    [r c] = size(senDatR);
    
    epi1DataRS = sum(senDatR(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
    epi1DataGS = sum(senDatG(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
    epi1DataBS = sum(senDatB(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);

    epi6DataRS = sum(senDatR(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
    epi6DataGS = sum(senDatG(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
    epi6DataBS = sum(senDatB(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);

    figure;
    hold on; 
    plot(epi1DataRS,'LineWidth',2,'color','red');
    plot(epi1DataGS,'LineWidth',2,'color','green');
    plot(epi1DataBS,'LineWidth',2,'color','blue');
    %plot(epi6DataRS,'LineWidth',2,'color','#A2142F');
    %plot(epi6DataGS,'LineWidth',2,'color','#77AC30');
    %plot(epi6DataBS,'LineWidth',2,'color','#0072BD');
    title('Epidermis Type I Sensor Averages')
    legend({'Red Pixel','Green Pixel','Blue Pixel'},'location','southwest')

    %{
    figure;
    hold on;
    plot(epi1DataGS,'LineWidth',2,'color','green');
    plot(epi6DataGS,'LineWidth',2,'color','#77AC30');
    title('Epidermis Type I Average Green Pixel Value');
    legend({'Large Blood Volume','Small Blood Volume'});
    %}

    figure;
    hold on;
    plot(epi6DataRS,'LineWidth',2,'color','red');
    plot(epi6DataGS,'LineWidth',2,'color','green');
    plot(epi6DataBS,'LineWidth',2,'color','blue');
    title('Epidermis Type VI Sensor Averages')
    legend({'Red Pixel','Green Pixel','Blue Pixel'},'location','southwest')
end

%bayer-1
if(1)
    senDat = sensor.data.volts;
    [r c] = size(senDat);
    senDatC = senDat(1:2:r,1:2:c);
    senDatM = senDat(2:2:r,2:2:c);
    senDatY = senDat(1:2:r,2:2:c) + senDat(2:2:r,1:2:c);
    figure; imagesc(senDatC); colorbar;
    figure; imagesc(senDatM); colorbar;
    figure; imagesc(senDatY); colorbar;
    [r c] = size(senDatC);
    
    epi1DataCS = sum(senDatC(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
    epi1DataYS = sum(senDatY(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
    epi1DataMS = sum(senDatM(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
   
    epi6DataCS = sum(senDatC(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
    epi6DataYS = sum(senDatY(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
    epi6DataMS = sum(senDatM(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
    

    figure;
    hold on; 
    plot(epi1DataCS,'LineWidth',2,'color','cyan');
    plot(epi1DataYS,'LineWidth',2,'color',"#EDB120");
    plot(epi1DataMS,'LineWidth',2,'color','magenta');
    title('Epidermis Type I Average Pixel Value');
    legend({'cyan','yellow','magenta'},'location','northeast')

    figure;
    hold on;
    plot(epi6DataCS,'LineWidth',2,'color','cyan');
    plot(epi6DataYS,'LineWidth',2,'color',"#EDB120");
    plot(epi6DataMS,'LineWidth',2,'color','magenta');
    title('Epidermis Type VI Average Pixel Value');
    legend({'cyan','yellow','magenta'},'location','northeast')

end

%Ar0132at RCCC
if(0)
    senDat = sensor.data.volts;
    [r c] = size(senDat);
    senDatR = senDat(2:2:r,2:2:c);
    senDatC = (senDat(1:2:r,2:2:c) + senDat(2:2:r,1:2:c)+senDat(1:2:r,1:2:c))/3;
    figure; imagesc(senDatC); colorbar;
    figure; imagesc(senDatR); colorbar;
    [r c] = size(senDatC);
    
    epi1DataCS = sum(senDatC(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
    epi1DataRS = sum(senDatR(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
   
    epi6DataCS = sum(senDatC(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
    epi6DataRS = sum(senDatR(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);
    

    figure;
    hold on; 
    plot(epi1DataCS,'LineWidth',2,'color','black');
    plot(epi1DataRS,'LineWidth',2,'color','red');
    title('Epidermis Type I Average Pixel Value');
    legend({'Clear Pixels', 'Red Pixels'},'Location','southwest')

    figure;
    hold on;
    plot(epi6DataCS,'LineWidth',2,'color','black');
    plot(epi6DataRS,'LineWidth',2,'color','red');
    title('Epidermis Type VI Average Pixel Value');
    legend({'Clear Pixels', 'Red Pixels'},'location','southwest')

end


%{
dataSensor = sensorGet(sensor,'volts');
epi1DataRS = sum(dataSensor(1+1:round(r/2)-1,:,1))/length(1+1:round(r/2)-1);
epi6DataRS = sum(dataSensor(round(r/2)+2:r-1,:,1))/length(round(r/2)+2:r-1);

figure(102);
hold on;
plot(epi1DataRS,'LineWidth',2);

figure(103);
hold on;
plot(epi6DataRS,'LineWidth',2);
%}