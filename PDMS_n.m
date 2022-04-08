function result=PDMS_n(x)
% 1) X. Zhang, J. Qiu, X. Li, J. Zhao, L. Liu. Complex refractive indices measurements of polymers in visible and near-infrared bands, Appl. Opt. 59, 2337-2344 (2020) (0.4-2 µm)
% 2) X. Zhang, J. Qiu, J. Zhao, X. Li, L. Liu. Complex refractive indices measurements of polymers in infrared bands, J. Quant. Spectrosc. Radiat. Transf. 252, 107063 (2020) (2-20 µm)
x=10^6*x;
data=[0.29999	1.42109
0.40	1.42109
0.41	1.41984
0.42	1.41859
0.43	1.41780
0.44	1.41642
0.45	1.41581
0.46	1.41499
0.47	1.41395
0.48	1.41342
0.49	1.41243
0.50	1.41189
0.51	1.41128
0.52	1.41051
0.53	1.41000
0.54	1.40954
0.55	1.40929
0.56	1.40851
0.57	1.40819
0.58	1.40776
0.59	1.40739
0.60	1.40715
0.61	1.40653
0.62	1.40620
0.63	1.40553
0.64	1.40527
0.65	1.40507
0.66	1.40484
0.67	1.40451
0.68	1.40453
0.69	1.40389
0.70	1.40397
0.71	1.40335
0.72	1.40324
0.73	1.40313
0.74	1.40307
0.75	1.40234
0.76	1.40235
0.77	1.40206
0.78	1.40201
0.79	1.40158
0.80	1.40124
0.81	1.40126
0.82	1.40150
0.83	1.40150
0.84	1.40137
0.85	1.40153
0.86	1.40137
0.87	1.40120
0.88	1.40091
0.89	1.40086
0.90	1.40077
0.91	1.40061
0.92	1.40047
0.93	1.40040
0.94	1.40025
0.95	1.40015
0.96	1.40015
0.97	1.39997
0.98	1.39977
0.99	1.39970
1.00	1.39976
1.01	1.39982
1.02	1.39993
1.03	1.39734
1.04	1.39708
1.05	1.39739
1.06	1.39685
1.07	1.39677
1.08	1.39696
1.09	1.39675
1.10	1.39673
1.11	1.39660
1.12	1.39654
1.13	1.39622
1.14	1.39568
1.15	1.39599
1.16	1.39631
1.17	1.39619
1.18	1.39656
1.19	1.39597
1.20	1.39567
1.21	1.39517
1.22	1.39569
1.23	1.39558
1.24	1.39494
1.25	1.39481
1.26	1.39515
1.27	1.39555
1.28	1.39490
1.29	1.39445
1.30	1.39483
1.31	1.39424
1.32	1.39390
1.33	1.39453
1.34	1.39340
1.35	1.39398
1.36	1.39455
1.37	1.39453
1.38	1.39425
1.39	1.39359
1.40	1.39387
1.41	1.39361
1.42	1.39392
1.43	1.39373
1.44	1.39366
1.45	1.39347
1.46	1.39318
1.47	1.39373
1.48	1.39400
1.49	1.39380
1.50	1.39423
1.51	1.39279
1.52	1.39283
1.53	1.39321
1.54	1.39284
1.55	1.39328
1.56	1.39309
1.57	1.39291
1.58	1.39231
1.59	1.39276
1.60	1.39251
1.61	1.39299
1.62	1.39257
1.63	1.39185
1.64	1.39271
1.65	1.39170
1.66	1.39185
1.67	1.39190
1.68	1.39262
1.69	1.39197
1.70	1.39116
1.71	1.39032
1.72	1.39027
1.73	1.39079
1.74	1.39191
1.75	1.39224
1.76	1.39094
1.77	1.39111
1.78	1.39197
1.79	1.39189
1.80	1.39084
1.81	1.39131
1.82	1.39026
1.83	1.39172
1.84	1.39211
1.85	1.39078
1.86	1.39137
1.87	1.39159
1.88	1.39150
1.89	1.39093
1.90	1.39098
1.91	1.39019
1.92	1.39004
1.93	1.39071
1.94	1.39044
1.95	1.39051
1.96	1.39081
1.97	1.39043
1.98	1.39098
1.99	1.38973
2.00	1.39092
2.0097	1.40028
2.0190	1.39960
2.0285	1.39855
2.0381	1.39720
2.0478	1.39612
2.0575	1.39699
2.0673	1.39903
2.0773	1.39881
2.0873	1.39687
2.0975	1.39696
2.1077	1.39914
2.1180	1.39967
2.1285	1.39987
2.1390	1.40233
2.1496	1.40370
2.1568	1.40326
2.1676	1.40133
2.1785	1.39864
2.1896	1.39930
2.1970	1.40095
2.2082	1.40247
2.2196	1.40440
2.2272	1.40536
2.2387	1.40374
2.2465	1.40253
2.2582	1.40081
2.2661	1.40030
2.2781	1.40186
2.2861	1.40307
2.2983	1.40251
2.3065	1.40157
2.3188	1.40147
2.3272	1.40085
2.3398	1.40004
2.3482	1.40030
2.3568	1.40021
2.3697	1.40048
2.3784	1.40061
2.3872	1.40076
2.3960	1.40089
2.4093	1.40133
2.4183	1.40200
2.4274	1.40179
2.4365	1.40081
2.4457	1.40033
2.4596	1.40147
2.4690	1.40204
2.4784	1.40209
2.4880	1.40105
2.4975	1.40027
2.5072	1.40027
2.5169	1.39981
2.5268	1.39979
2.5366	1.40073
2.5466	1.40055
2.5567	1.39877
2.5668	1.39652
2.5770	1.39631
2.5873	1.39760
2.5976	1.39832
2.6081	1.39925
2.6186	1.40041
2.6293	1.39945
2.6400	1.39704
2.6454	1.39654
2.6562	1.39786
2.6671	1.39947
2.6782	1.39874
2.6893	1.39730
2.6949	1.39750
2.7061	1.39903
2.7175	1.39970
2.7289	1.39986
2.7347	1.40003
2.7462	1.40070
2.7579	1.39972
2.7697	1.39755
2.7756	1.39727
2.7876	1.39815
2.7996	1.39933
2.8057	1.39987
2.8179	1.39966
2.8240	1.39890
2.8364	1.39710
2.8489	1.39627
2.8551	1.39635
2.8678	1.39733
2.8741	1.39768
2.8869	1.39797
2.8998	1.39849
2.9063	1.39860
2.9194	1.39795
2.9260	1.39761
2.9393	1.39743
2.9460	1.39741
2.9594	1.39716
2.9662	1.39669
2.9798	1.39594
2.9867	1.39628
2.9936	1.39706
3.0075	1.39749
3.0145	1.39719
3.0286	1.39643
3.0357	1.39647
3.0499	1.39718
3.0571	1.39727
3.0644	1.39740
3.0789	1.39742
3.0863	1.39713
3.0936	1.39680
3.1085	1.39584
3.1159	1.39496
3.1234	1.39379
3.1386	1.39324
3.1462	1.39379
3.1538	1.39429
3.1693	1.39476
3.1770	1.39433
3.1848	1.39372
3.1927	1.39301
3.2085	1.39247
3.2164	1.39263
3.2244	1.39267
3.2325	1.39253
3.2487	1.39146
3.2569	1.39121
3.2651	1.39100
3.2733	1.39072
3.2899	1.39062
3.2983	1.38988
3.3067	1.38851
3.3152	1.38681
3.3237	1.38496
3.3322	1.38255
3.3494	1.37358
3.3581	1.37134
3.3668	1.37853
3.3756	1.39535
3.3844	1.41081
3.3933	1.41582
3.4022	1.41294
3.4111	1.40856
3.4292	1.40376
3.4383	1.40340
3.4474	1.40403
3.4566	1.40436
3.4658	1.40409
3.4751	1.40326
3.4845	1.40221
3.4939	1.40125
3.5033	1.40063
3.5128	1.40039
3.5224	1.40016
3.5320	1.39951
3.5416	1.39853
3.5513	1.39792
3.5611	1.39841
3.5709	1.39924
3.5807	1.39928
3.5907	1.39851
3.6006	1.39810
3.6107	1.39884
3.6207	1.39913
3.6309	1.39830
3.6411	1.39753
3.6513	1.39735
3.6617	1.39690
3.6720	1.39586
3.6825	1.39530
3.6930	1.39548
3.7035	1.39588
3.7141	1.39610
3.7248	1.39566
3.7355	1.39495
3.7463	1.39450
3.7572	1.39453
3.7681	1.39495
3.7791	1.39499
3.7901	1.39450
3.8013	1.39397
3.8124	1.39377
3.8237	1.39349
3.8350	1.39339
3.8464	1.39371
3.8578	1.39355
3.8693	1.39276
3.8809	1.39219
3.8926	1.39184
3.9043	1.39149
3.9161	1.39113
3.9280	1.39115
3.9399	1.39175
3.9519	1.39239
3.9640	1.39221
3.9762	1.39149
3.9884	1.39120
4.0007	1.39128
4.0131	1.39093
4.0255	1.39025
4.0381	1.38993
4.0507	1.38997
4.0634	1.39025
4.0762	1.39022
4.0890	1.38974
4.1020	1.38930
4.1150	1.38929
4.1281	1.38936
4.1413	1.38912
4.1546	1.38879
4.1679	1.38875
4.1814	1.38874
4.1949	1.38886
4.2085	1.38916
4.2222	1.38994
4.2360	1.38929
4.2499	1.38714
4.2639	1.38749
4.2780	1.38893
4.2921	1.38864
4.3064	1.38701
4.3208	1.38723
4.3352	1.38775
4.3498	1.38746
4.3644	1.38661
4.3791	1.38593
4.3940	1.38561
4.4089	1.38551
4.4240	1.38542
4.4391	1.38535
4.4544	1.38523
4.4697	1.38444
4.4852	1.38377
4.5008	1.38368
4.5165	1.38353
4.5323	1.38336
4.5482	1.38332
4.5642	1.38294
4.5803	1.38233
4.5965	1.38243
4.6129	1.38292
4.6294	1.38280
4.6460	1.38261
4.6627	1.38250
4.6795	1.38226
4.6965	1.38173
4.7136	1.38121
4.7308	1.38039
4.7481	1.37976
4.7655	1.37976
4.7831	1.37999
4.8008	1.37979
4.8187	1.37925
4.8367	1.37867
4.8548	1.37801
4.8730	1.37747
4.8914	1.37730
4.9100	1.37710
4.9286	1.37647
4.9474	1.37595
4.9664	1.37552
4.9855	1.37505
5.0047	1.37520
5.0241	1.37538
5.0437	1.37499
5.0634	1.37463
5.0832	1.37450
5.1033	1.37445
5.1234	1.37467
5.1438	1.37432
5.1643	1.37384
5.1849	1.37399
5.2057	1.37398
5.2267	1.37307
5.2479	1.37199
5.2692	1.37158
5.2907	1.37117
5.3124	1.37023
5.3343	1.36935
5.3563	1.36952
5.3785	1.36938
5.4009	1.36887
5.4235	1.36838
5.4463	1.36772
5.4693	1.36709
5.4925	1.36624
5.5159	1.36549
5.5394	1.36500
5.5632	1.36410
5.5872	1.36347
5.6114	1.36311
5.6358	1.36284
5.6604	1.36302
5.6852	1.36229
5.7103	1.36127
5.7355	1.36110
5.7610	1.36045
5.7867	1.35957
5.8127	1.35862
5.8389	1.35789
5.8653	1.35815
5.8919	1.35778
5.9188	1.35573
5.9460	1.35404
5.9734	1.35385
6.0011	1.35358
6.0290	1.35284
6.0571	1.35203
6.0856	1.35100
6.1143	1.35000
6.1433	1.34880
6.1725	1.34723
6.2020	1.34606
6.2319	1.34605
6.2620	1.34595
6.2924	1.34520
6.3231	1.34372
6.3541	1.34273
6.3854	1.34234
6.4170	1.34110
6.4489	1.33891
6.4811	1.33700
6.5137	1.33512
6.5466	1.33371
6.5798	1.33250
6.6134	1.33083
6.6473	1.32824
6.6816	1.32662
6.7162	1.32489
6.7512	1.32235
6.7865	1.31995
6.8223	1.31759
6.8583	1.31451
6.8948	1.31220
6.9317	1.31140
6.9690	1.30858
7.0066	1.30449
7.0447	1.30322
7.0832	1.30571
7.1221	1.30963
7.1615	1.31189
7.2013	1.30840
7.2415	1.30095
7.2822	1.29452
7.3233	1.28974
7.3649	1.28464
7.4070	1.27893
7.4496	1.27215
7.4926	1.26451
7.5362	1.25588
7.5803	1.24676
7.6249	1.23572
7.6700	1.22274
7.7156	1.20633
7.7618	1.18428
7.8086	1.13732
7.8559	1.06690
7.9038	1.13149
7.9523	1.31562
8.0014	1.42755
8.0511	1.37931
8.1014	1.31645
8.1524	1.27762
8.2040	1.24968
8.2562	1.22830
8.3091	1.21345
8.3628	1.20275
8.4171	1.19147
8.4721	1.17666
8.5278	1.15959
8.5843	1.14061
8.6415	1.11362
8.6995	1.07813
8.7583	1.03538
8.8179	0.98636
8.8783	0.93665
8.9395	0.90189
9.0016	0.91161
9.0645	0.99110
9.1284	1.12930
9.1931	1.28201
9.2588	1.40749
9.3254	1.49129
9.3930	1.54005
9.4615	1.56597
9.5311	1.57798
9.6017	1.58843
9.6733	1.61404
9.7461	1.67584
9.8199	1.77617
9.8949	1.85190
9.9710	1.84298
10.048	1.78062
10.127	1.71731
10.207	1.66568
10.288	1.62340
10.370	1.58894
10.453	1.55986
10.538	1.53370
10.625	1.50932
10.713	1.48704
10.802	1.46405
10.893	1.43983
10.985	1.41698
11.079	1.39386
11.174	1.36267
11.272	1.32503
11.370	1.29663
11.471	1.30210
11.573	1.32967
11.678	1.33658
11.784	1.33938
11.892	1.34524
12.002	1.29807
12.114	1.24035
12.229	1.26473
12.345	1.38610
12.464	1.61522
12.585	1.86651
12.708	1.96260
12.834	1.89309
12.962	1.79395
13.093	1.73241
13.227	1.71030
13.363	1.69424
13.502	1.66425
13.645	1.63582
13.790	1.61199
13.938	1.59044
14.089	1.57874
14.244	1.57910
14.403	1.58302
14.564	1.58508
14.730	1.58381
14.899	1.57887
15.072	1.57711
15.250	1.57350
15.431	1.56094
15.617	1.54848
15.808	1.54164
16.003	1.53726
16.203	1.52800
16.408	1.51659
16.618	1.51254
16.834	1.51107
17.056	1.50150
17.283	1.49402
17.517	1.49619
17.757	1.50036
18.003	1.49790
18.257	1.49605
18.518	1.49709
18.786	1.49362
19.062	1.48785
19.347	1.48798
19.640	1.48985
19.942	1.49006
30.942	1.49006];
result=interp1(data(:,1),data(:,2),x);