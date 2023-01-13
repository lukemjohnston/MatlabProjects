%uiimport('spectrum.xls');
x=spectrum{:,1};
y=spectrum{:,2};
plot(x,y);


[s1, s2, s3, s4, s5, s6] = deal(4000);   %start variables
[f1, f2, f3, f4, f5, f6] = deal(0);      %finish variables

for i=1:4000
    if ((x(i) > 3.669e+14) && (x(i) < 3.673e+14))       %line #1
        if (i<s1)
            s1=i;
        end
        if (i>f1)
            f1=i;
        end
    end

    if ((x(i) > 3.774e+14) && (x(i) < 3.778e+14))       %line #2
        if (i<s2)
            s2=i;
        end
        if (i>f2)
            f2=i;
        end
    end

    if ((x(i) > 3.82e+14) && (x(i) < 3.824e+14))       %line #3
        if (i<s3)
            s3=i;
        end
        if (i>f3)
            f3=i;
        end
    end

    if ((x(i) > 3.961e+14) && (x(i) < 3.965e+14))       %line #4
        if (i<s4)
            s4=i;
        end
        if (i>f4)
            f4=i;
        end
    end

    if ((x(i) > 3.9686e+14) && (x(i) < 3.9726e+14))       %line #5
        if (i<s5)
            s5=i;
        end
        if (i>f5)
            f5=i;
        end
    end

    if ((x(i) > 4.224e+14) && (x(i) < 4.228e+14))       %line #6
        if (i<s6)
            s6=i;
        end
        if (i>f6)
            f6=i;
        end
    end
end


Qsum1 = 0;          %calculate the sum of trapezoids for the spectral line #1
for i=s1:f1-1   
    Q = 1/2*(x(i+1)-x(i))*(y(i)+y(i+1));
    Qsum1 = Qsum1 + Q;
end
%plot(x(s1:f1),y(s1:f1));
disp('strength of spectral line #1:');
disp(Qsum1);


Qsum2 = 0;          %calculate the sum of trapezoids for the spectral line #2
for i=s2:f2-1
    Q = 1/2*(x(i+1)-x(i))*(y(i)+y(i+1));
    Qsum2 = Qsum2 + Q;
end
%plot(x(s2:f2),y(s2:f2));
disp('strength of spectral line #2:');
disp(Qsum2);


Qsum3 = 0;          %calculate the sum of trapezoids for the spectral line #3
for i=s3:f3-1
    Q = 1/2*(x(i+1)-x(i))*(y(i)+y(i+1));
    Qsum3 = Qsum3 + Q;
end
%plot(x(s3:f3),y(s3:f3));
disp('strength of spectral line #3:');
disp(Qsum3);


Qsum4 = 0;          %calculate the sum of trapezoids for the spectral line #4
for i=s4:f4-1
    Q = 1/2*(x(i+1)-x(i))*(y(i)+y(i+1));
    Qsum4 = Qsum4 + Q;
end
%plot(x(s4:f4),y(s4:f4));
disp('strength of spectral line #4:');
disp(Qsum4);


Qsum5 = 0;          %calculate the sum of trapezoids for the spectral line #5
for i=s5:f5-1
    Q = 1/2*(x(i+1)-x(i))*(y(i)+y(i+1));
    Qsum5 = Qsum5 + Q;
end
%plot(x(s5:f5),y(s5:f5));
disp('strength of spectral line #5:');
disp(Qsum5);


Qsum6 = 0;          %calculate the sum of trapezoids for the spectral line #5
for i=s6:f6-1
    Q = 1/2*(x(i+1)-x(i))*(y(i)+y(i+1));
    Qsum6 = Qsum6 + Q6;
end
%plot(x(s6:f6),y(s6:f6));
disp('strength of spectral line #6:');
disp(Qsum6);


title('Spectral Lines')
xlabel('Frequency Hz')
ylabel('Specific Intensity Hz')

