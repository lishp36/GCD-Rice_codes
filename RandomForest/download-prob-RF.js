
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
    yr, dStart2, dEnd2,
    studyPlace, pre, dStep2, unit2, mask, band2
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
      // .map(pre.renameBandsS2)
      .map(function temp (img) {return img.select('B11').rename(band2)})
      //.map(pre.getNDVI).map(pre.getLSWI)//.map(pre.getREP).map(pre.getGCVI)
      ;
    var collection_L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LT05/C02/T2_L2"))
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
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL8C2)
      .map(pre.rmCloudLC2).map(pre.rmCloudShadowLC2).map(pre.rmSnowLC2)
      .map(pre.getNDVI)//.map(pre.getGCVI)
      ;
    var collection_L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LC09/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL8C2)
      .map(pre.rmCloudLC2).map(pre.rmCloudShadowLC2).map(pre.rmSnowLC2)
      .map(pre.getNDVI)//.map(pre.getGCVI)
      ;
    var collection_opt =
      collection_S2
      // collection_L8
      .merge(collection_L5)
      .merge(collection_L7)
      .merge(collection_L8)
      .merge(collection_L9)
      ;
  
    // composte
    var size2 = dEnd2.difference(dStart2, unit2).divide(dStep2).toInt().getInfo();
    var list2 = ee.List.sequence(0, size2-1, 1);
    var collection2 = ee.ImageCollection(list2.map(function (d) {
      var d_start = dStart2.advance(ee.Number(d).multiply(dStep2), unit2);
      var d_end = d_start.advance(dStep2, unit2);
      var collection_opt_ = collection_opt.select(band2).filterDate(d_start, d_end)
        .select(band2).median().clip(studyPlace);
      return ee.Image(ee.Algorithms.If({
        condition: collection_opt_.bandNames().length().neq(1),
        trueCase: ee.Image.constant(0).rename(band2),
        falseCase: collection_opt_
      }));
    }));
  
    var collection2_ = collection2;
  
    var result = collection2_.toBands();
    return result.clip(studyPlace).unmask().updateMask(mask);
  };
  
  var pre = require("path/to/work:preprocess");
  var At = "place"; var province = At;
  var chinaProvinces = ee.FeatureCollection("path/to/shapefile");
  var studyPlace = chinaProvinces.filterMetadata("", "equals", "");
  
  var version = "-vrf7";
  var select = ee.FeatureCollection(
    "path/" + At + "-S2_L8_L9_swir1-1_8_361_day"
  );
  
  var band2 = ["swir1"];
  Map.addLayer(studyPlace, {}, "study place");
  
  for (var yr = 1990; yr <= 2023; yr++) {
    var start2 = ""; var end2 = ""; var dStep2 = 8; var unit2 = "day";
    var dStart2 = ee.Date.parse("yyyyDDD", yr + start2);
    var dEnd2 = ee.Date.parse("yyyyDDD", yr + end2).advance(dStep2, unit2);
    
    var mask = ee.Image.constant(1);
    
    var compositedImageClip = exports.getCollection(
      yr, dStart2, dEnd2,
      studyPlace, pre, dStep2, unit2, mask, band2
    ).subtract(0).multiply(1000).toInt16()//.rename(bands);
    
    print("compositedImageClip", compositedImageClip);
    // Map.addLayer(compositedImageClip.clip(studyPlace), {}, "Image");
  
    var trainargs = {
      features: select,
      classProperty: "class",
      inputProperties: compositedImageClip.bandNames()
    };
    
    var trained = ee.Classifier.smileRandomForest({
      numberOfTrees: 100,
      seed: 9999
    }).setOutputMode("MULTIPROBABILITY").train(trainargs);
    
    var classified = ee.Image(
      compositedImageClip.classify(trained)
    ).arrayGet(1).clip(studyPlace);
    
    print(classified);
    
    Export.image.toDrive({
      image: classified,
      description: "prob-" + At + "-" + yr + "-rice" + version,
      scale: 30,
      fileNamePrefix: "prob-" + At + "-" + yr + "-rice" + version + "-WGS84",
      maxPixels: 1e13,
      region: studyPlace,
      folder: "",
      crs: "EPSG:4326",
    });
  }
  
  