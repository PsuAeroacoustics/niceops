stk.v.12.0
WrittenBy    STK_v12.0.0

BEGIN ReportStyle

    BEGIN ClassId
        Class		 Aircraft
    END ClassId

    BEGIN Header
        StyleType		 0
        Title		 LLA Position & Velocity
        Date		 Yes
        Name		 Yes
        IsHidden		 No
        DescShort		 No
        DescLong		 No
        YLog10		 No
        Y2Log10		 No
        YUseWholeNumbers		 No
        Y2UseWholeNumbers		 No
        VerticalGridLines		 No
        HorizontalGridLines		 No
        AnnotationType		 Spaced
        NumAnnotations		 3
        NumAngularAnnotations		 5
        ShowYAnnotations		 Yes
        AnnotationRotation		 1
        BackgroundColor		 #ffffff
        ForegroundColor		 #000000
        ViewableDuration		 0
        RealTimeMode		 No
        DayLinesStatus		 1
        LegendStatus		 1
        LegendLocation		 1

        BEGIN PostProcessor
            Destination		 0
            Use		 0
            Destination		 1
            Use		 0
            Destination		 2
            Use		 0
            Destination		 3
            Use		 0
        END PostProcessor
        NumSections		 1
    END Header

    BEGIN Section
        Name		 Section 1
        ClassName		 Aircraft
        NameInTitle		 No
        ExpandMethod		 0
        PropMask		 2
        ShowIntervals		 No
        NumIntervals		 0
        NumLines		 1

        BEGIN Line
            Name		 Line 1
            NumElements		 17

            BEGIN Element
                Name		 Time
                IsIndepVar		 Yes
                IndepVarName		 Time
                Title		 Time
                NameInTitle		 No
                Service		 LLAState
                Type		 Fixed
                Element		 Time
                SumAllowedMask		 0
                SummaryOnly		 No
                DataType		 0
                UnitType		 2
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 LLA State-Fixed-Lat
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Lat
                NameInTitle		 No
                Service		 LLAState
                Type		 Fixed
                Element		 Lat
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 19
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 No
                BEGIN Units
                    LatitudeUnit		 Degrees
                END Units
            END Element

            BEGIN Element
                Name		 LLA State-Fixed-Lon
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Lon
                NameInTitle		 No
                Service		 LLAState
                Type		 Fixed
                Element		 Lon
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 20
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 No
                BEGIN Units
                    LongitudeUnit		 Degrees
                END Units
            END Element

            BEGIN Element
                Name		 LLA State-Fixed-Alt
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Alt
                NameInTitle		 No
                Service		 LLAState
                Type		 Fixed
                Element		 Alt
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 0
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 LLA State-Fixed-Lat Rate
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Lat Rate
                NameInTitle		 No
                Service		 LLAState
                Type		 Fixed
                Element		 Lat Rate
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 7
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 LLA State-Fixed-Lon Rate
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Lon Rate
                NameInTitle		 No
                Service		 LLAState
                Type		 Fixed
                Element		 Lon Rate
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 7
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 LLA State-Fixed-Alt Rate
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Alt Rate
                NameInTitle		 No
                Service		 LLAState
                Type		 Fixed
                Element		 Alt Rate
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 4
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 Heading-TrueOfDateRotating-Velocity
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Velocity
                NameInTitle		 Yes
                Service		 Heading
                Type		 TrueOfDateRotating
                Element		 Velocity
                SumAllowedMask		 1543
                SummaryOnly		 No
                DataType		 0
                UnitType		 4
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 Flight Attitude-Yaw
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Yaw
                NameInTitle		 Yes
                Service		 FlightAttitude
                Element		 Yaw
                SumAllowedMask		 1543
                SummaryOnly		 No
                DataType		 0
                UnitType		 20
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 Flight Attitude-Pitch
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Pitch
                NameInTitle		 Yes
                Service		 FlightAttitude
                Element		 Pitch
                SumAllowedMask		 1543
                SummaryOnly		 No
                DataType		 0
                UnitType		 19
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 Flight Attitude-Roll
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Roll
                NameInTitle		 Yes
                Service		 FlightAttitude
                Element		 Roll
                SumAllowedMask		 1543
                SummaryOnly		 No
                DataType		 0
                UnitType		 3
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 Flight Attitude-Bank
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Bank
                NameInTitle		 Yes
                Service		 FlightAttitude
                Element		 Bank
                SumAllowedMask		 1543
                SummaryOnly		 No
                DataType		 0
                UnitType		 3
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 Flight Attitude-Heading
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 Heading
                NameInTitle		 Yes
                Service		 FlightAttitude
                Element		 Heading
                SumAllowedMask		 1543
                SummaryOnly		 No
                DataType		 0
                UnitType		 3
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 Yes
            END Element

            BEGIN Element
                Name		 Cartesian Acceleration-Fixed-x
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 ax
                NameInTitle		 Yes
                Service		 CartAcc
                Type		 Fixed
                Element		 x
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 11
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 No
                BEGIN Units
                    DistanceUnit		 Meters
                    TimeUnit		 Seconds
                END Units
            END Element

            BEGIN Element
                Name		 Cartesian Acceleration-Fixed-y
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 ay
                NameInTitle		 Yes
                Service		 CartAcc
                Type		 Fixed
                Element		 y
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 11
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 No
                BEGIN Units
                    DistanceUnit		 Meters
                    TimeUnit		 Seconds
                END Units
            END Element

            BEGIN Element
                Name		 Cartesian Acceleration-Fixed-z
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 az
                NameInTitle		 Yes
                Service		 CartAcc
                Type		 Fixed
                Element		 z
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 11
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 No
                BEGIN Units
                    DistanceUnit		 Meters
                    TimeUnit		 Seconds
                END Units
            END Element

            BEGIN Element
                Name		 Cartesian Acceleration-Fixed-magnitude
                IsIndepVar		 No
                IndepVarName		 Time
                Title		 acc mag
                NameInTitle		 Yes
                Service		 CartAcc
                Type		 Fixed
                Element		 magnitude
                SumAllowedMask		 1559
                SummaryOnly		 No
                DataType		 0
                UnitType		 11
                LineStyle		 0
                LineWidth		 0
                PointStyle		 0
                PointSize		 0
                FillPattern		 0
                LineColor		 #000000
                FillColor		 #000000
                PropMask		 0
                UseScenUnits		 No
                BEGIN Units
                    DistanceUnit		 Meters
                    TimeUnit		 Seconds
                END Units
            END Element
        END Line
    END Section

    BEGIN LineAnnotations
    END LineAnnotations
END ReportStyle

