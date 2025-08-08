# Sample Data Directory

Bu klasör, deprem sinyali analiz GUI'sini test etmek için kullanılacak örnek veri dosyalarını içerir.

## Dosya Formatı

GUI uygulaması MATLAB (.mat) formatındaki dosyaları bekler. Her dosya aşağıdaki yapıya sahip olmalıdır:

```matlab
EQ.anEQ.Accel      % 3 kanallı ivme verisi (North, East, Up)
EQ.anEQ.Ptime      % P dalgası geliş zamanı
EQ.anEQ.Stime      % S dalgası geliş zamanı
EQ.anEQ.epicenter  % Epicenter koordinatları [enlem, boylam]
EQ.anEQ.statco     % İstasyon koordinatları [enlem, boylam]
EQ.anEQ.place      % Deprem yeri
EQ.anEQ.date       % Deprem tarihi
EQ.anEQ.depth      % Derinlik (km)
EQ.anEQ.magnitudeval % Büyüklük
EQ.anEQ.magnitudetype % Büyüklük tipi
EQ.anEQ.statID     % İstasyon ID
EQ.anEQ.numofData  % Veri noktası sayısı
EQ.anEQ.pga        % Peak Ground Acceleration
EQ.anEQ.Vs30       % Vs30 değeri
```

## Örnek Dosyalar

Bu klasöre test etmek istediğiniz .mat dosyalarını ekleyebilirsiniz. Dosyalar otomatik olarak GUI'de listelenecektir.

## Not

- Dosyalar .mat formatında olmalıdır
- Her dosya yukarıdaki yapıya uygun olmalıdır
- GUI uygulaması bu klasördeki tüm .mat dosyalarını otomatik olarak yükleyecektir
