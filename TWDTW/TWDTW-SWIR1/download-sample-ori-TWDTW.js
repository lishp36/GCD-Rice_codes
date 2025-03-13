// download-sample-ori-twdtw.js
// 下载样本点上的时间序列

/**
 * 指定时间段和波段，获取 SAR 和 光学（S2, L7, L8）影像
 * @param dStart1 SAR 开始日期
 * @param dEnd1 SAR 结束日期
 * @param dStart2 光学开始日期
 * @param dEnd2 光学结束日期
 * @param {ee.Feature} studyPlace 研究区域的范围
 * @param {Function} pre 预处理函数
 * @param dStep1 SAR 合成（日/月数）
 * @param unit1 SAR 合成单位 day/month
 * @param dStep2 光学合成（日/月数）
 * @param unit2 光学合成单位 day/month
 * @param {ee.Image} mask 耕地掩膜
 * @param {List} band1 SAR 波段
 * @param {List} band2 光学波段
 * @returns 合成的影像（时间合成到波段）
 */
exports.getCollection = function (
    yr, dStart1, dEnd1, dStart2, dEnd2,
    studyPlace, pre, dStep1, unit1, dStep2, unit2, mask, band1, band2, description
  ) {
    var dateFilter2 = ee.Filter.date(
      dStart2.advance(-dStep2, unit2), dEnd2.advance(dStep2, unit2)
    );
    var S2 = ee.ImageCollection('COPERNICUS/S2')
      .filter(dateFilter2).filterBounds(studyPlace);
    var S2C = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
      .filter(dateFilter2).filterBounds(studyPlace);
    var collection_S2 = pre.indexJoin(S2, S2C, 'cloud_probability')
      .map(pre.scaleFactorsS2).map(pre.maskImage_S2)//.map(pre.rmCloudS2)
      .map(pre.renameBandsS2)
      .map(pre.getNDVI).map(pre.getLSWI)//.map(pre.getREP).map(pre.getGCVI)
      ;
    var collection_L5 = ee.ImageCollection("LANDSAT/LE05/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LE05/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL7C2)
      .map(pre.rmCloudLC2).map(pre.rmCloudShadowLC2).map(pre.rmSnowLC2)
      .map(pre.getNDVI)//.map(pre.getGCVI)
      ;
    var collection_L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LE07/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL7C2)
      .map(pre.rmCloudLC2).map(pre.rmCloudShadowLC2).map(pre.rmSnowLC2)
      .map(pre.getNDVI)//.map(pre.getGCVI)
      ;
    var collection_L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LC08/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL8C2).map(pre.rmCloudLC2)
      .map(pre.getNDVI).map(pre.getGCVI)
      ;
    var collection_L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LC09/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL8C2)
      .map(pre.rmCloudLC2).map(pre.rmCloudShadowLC2).map(pre.rmSnowLC2)
      .map(pre.getNDVI).map(pre.getGCVI)
      ;
    var collection_opt =
      collection_S2
      // collection_L8
      //collection_L7
      // .merge(collection_L7)
      .merge(collection_L8)
      .merge(collection_L9)
      .select(band2)
      ;
  
    var dateFilter1 = ee.Filter.date(
      dStart1.advance(-dStep1, unit1), dEnd1.advance(dStep1, unit1)
    );
    var collection_S1 = ee.ImageCollection("COPERNICUS/S1_GRD")
      // .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
      .filter(ee.Filter.eq('instrumentMode', 'IW'))
      // .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
      .filter(dateFilter1).filterBounds(studyPlace)
      ;
  
    // composte
    var size2 = dEnd2.difference(dStart2, unit2).divide(dStep2).toInt().getInfo();
    var list2 = ee.List.sequence(0, size2-1, 1);
    var collection2 = ee.ImageCollection(list2.map(function (d) {
      var d_start = dStart2.advance(ee.Number(d).multiply(dStep2), unit2);
      var d_end = d_start.advance(dStep2, unit2);
      var collection_opt_ = collection_opt.filterDate(d_start, d_end)
        .select(band2).median().clip(studyPlace);
      return collection_opt_;
    }));
  
    var size1 = dEnd1.difference(dStart1, unit1).divide(dStep1).toInt().getInfo();
    var list1 = ee.List.sequence(0, size1-1, 1);
  
    var collectionS1 = ee.ImageCollection(list1.map(function (d) {
      var d_start = dStart1.advance(ee.Number(d).multiply(dStep1), unit1);
      var d_end = d_start.advance(dStep1, unit1);
      var collection_S1_ = collection_S1.filterDate(d_start, d_end)
        .median().select(band1).clip(studyPlace);
      return collection_S1_;
    }));
  
    var collection2_ = collection2;
    var collectionS1_ = collectionS1;
  
    if (band1.length > 0) {
      var result = collectionS1_.toBands().addBands(collection2_.toBands());
      return result.clip(studyPlace).unmask().updateMask(mask);
    } else {
      var result = collection2_.toBands().addBands(collectionS1_.toBands());
      return result.clip(studyPlace).unmask().updateMask(mask);
    }
  
  };
  
  var pre = require("");
  var cover = "rice"; // double_rice single_rice other_crop rice
    var At = ""; var province = At;
  
  var band1 = [
    // "VV",
    "VH",
  ];
  var band2 = [
    // "blue",
    // "green",
    // "red",
    // "nir",
    // "re1",
    // "re2",
    // "re3",
    // "re4",
    // "REP",
    // "swir1",
    // "swir2",
    // "NDVI",
    // "LSWI",
    // "GCVI",
  ];

  var points = ""
  var select = ee.FeatureCollection(points);  
  Map.addLayer(select, {}, "select");
  for (var yr = 2016; yr <= 2016; yr++) {
    var select1 = select;
    var chinaProvinces = ee.FeatureCollection("");
    var studyPlace = chinaProvinces.filterMetadata("", "equals", "");

   
    Map.addLayer(studyPlace, {}, "Image1");
     
    var mask = ee.Image.constant(1);
  
    var start1 = "1"; var end1 = "361"; var dStep1 = 12; var unit1 = "day";
    var dStart1 = ee.Date.parse("yyyyDDD", yr + start1);
    var dEnd1 = ee.Date.parse("yyyyDDD", yr + end1).advance(dStep1, unit1);
  
    var start2 = "1"; var end2 = "361"; var dStep2 = 8; var unit2 = "day";
    var dStart2 = ee.Date.parse("yyyyDDD", yr + start2);
    var dEnd2 = ee.Date.parse("yyyyDDD", yr + end2).advance(dStep2, unit2);
  
  
    var size1 = dEnd1.difference(dStart1, unit1).divide(dStep1).toInt().getInfo();
    var n = 0;
    var bands1 = [];
    for (var j = 0; j < band1.length; j++) {
      for (var i = 0; i < size1; i++) {
        bands1[n] = band1[j] + "_" + i.toString(); n += 1;
      }
    }
    var size2 = dEnd2.difference(dStart2, unit2).divide(dStep2).toInt().getInfo();
    var n = 0;
    var bands2 = [];
    for (var j = 0; j < band2.length; j++) {
      for (var i = 0; i < size2; i++) {
        bands2[n] = band2[j] + "_" + i.toString(); n += 1;
      }
    }
    var bands = bands1.concat(bands2);
    print(bands);
    var description = At + "-" + yr;
    if (band1.length > 0) {
      description = description + "-S1";
      for (var j = 0; j < band1.length; j++) {
        description = description + "_" + band1[j];
      }
      description = description + "-" + start1 + "_" + dStep1 + "_" + end1 + "_" + unit1;
    }
    if (band2.length > 0) {
      description = description + "-S2_L8_L9";
      for (var j = 0; j < band2.length; j++) {
        description = description + "_" + band2[j];
      }
      description = description + "-" + start2 + "_" + dStep2 + "_" + end2 + "_" + unit2;
    }
  
    var compositedImageClip = exports.getCollection(
      yr, dStart1, dEnd1, dStart2, dEnd2,
      studyPlace, pre, dStep1, unit1, dStep2, unit2, mask, band1, band2, description
    ).clip(studyPlace).subtract(0).multiply(100).toInt16().unmask();//.rename(bands)
  
    // print("compositedImageClip", compositedImageClip);
    Map.addLayer(compositedImageClip, {}, "Image");
  

    var select2 = compositedImageClip.sampleRegions({
      collection: select1, // maize: 1
      properties: ["class","b1","geometry"],//, "year"
      scale: 30,
      geometries: true,
    })
    ;
    print(select2);
    Export.table.toDrive({
      collection: select2,
      description: description + "_points_" + pre.dateTimeStr(),
      folder: "",
      fileNamePrefix: points
      //fileNamePrefix: description + "-" + cover
    });
  }
  