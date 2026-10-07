# İKARUS · ASRW 2026 Rehber ve Uçuş Simülatörü

Airbus Sloshing Rocket Workshop 2026 için Türkçe konu anlatımı ve su roketi uçuş simülatörü.
Site tamamen statiktir: sunucu tarafı kod, veritabanı veya derleme adımı gerekmez.

## Klasör yapısı

```
IKARUS-ASRW-Site/
├── index.html              Sayfanın tamamı (konu anlatımı + simülasyon arayüzü)
├── assets/
│   ├── css/style.css       Bütün stiller
│   ├── js/app.js           Fizik modeli, 2B/3B çizim, grafikler, rapor
│   └── img/favicon.svg     Sekme simgesi
├── robots.txt
├── .nojekyll               GitHub Pages için (dosyaları olduğu gibi yayınlar)
└── README.md
```

## Yayına alma

Klasörün **içeriğini** (index.html'in en üstte olduğu hâliyle) aşağıdakilerden birine yükleyin.

**GitHub Pages**
1. GitHub'da yeni bir depo açın ve bu klasördeki dosyaları deponun köküne yükleyin.
2. Settings → Pages → Source: "Deploy from a branch", Branch: `main`, klasör: `/ (root)`.
3. Birkaç dakika sonra site `https://KULLANICI-ADI.github.io/DEPO-ADI/` adresinde açılır.

**Netlify Drop** (hesap açıp sürükle-bırak)
1. https://app.netlify.com/drop adresine gidin.
2. Bu klasörü sayfaya sürükleyin; size bir adres verilir. Site adını Netlify ayarlarından değiştirebilirsiniz.

**Cloudflare Pages / Vercel**
"Static site" olarak yeni proje oluşturun, bu klasörü yükleyin. Build komutu yok, çıktı klasörü `/`.

**Kendi sunucunuz**
Dosyaları web sunucusunun kök klasörüne kopyalamanız yeterli (Apache, Nginx, IIS).

## Yerelde deneme

`index.html` dosyasına çift tıklayarak açabilirsiniz. Tarayıcıların dosya güvenlik kuralları
nedeniyle en doğrusu küçük bir yerel sunucu kullanmaktır:

```
cd IKARUS-ASRW-Site
python -m http.server 8000
```

Ardından tarayıcıda `http://localhost:8000` adresini açın.

## Dış bağlantılar

Sayfa şu kaynakları internetten yükler; bu yüzden ziyaretçinin internete bağlı olması gerekir:

- Chart.js 4.4.1 (grafikler) — cdn.jsdelivr.net
- three.js r128 (3B görünüm, yalnızca 3B seçilince) — cdn.jsdelivr.net, olmazsa cdnjs.cloudflare.com
- Yazı tipleri (Oxanium, IBM Plex Sans/Mono) — Google Fonts

Tarayıcının yerel depolaması (localStorage) yalnızca 2B/3B ve kamera tercihini hatırlamak için kullanılır.
Kişisel veri toplanmaz, sunucuya veri gönderilmez.

## Düzenleme

- Konu anlatımının metinleri `index.html` içinde, `<article class="topic" id="t-...">` blokları hâlindedir.
- Simülasyonun fizik modeli `assets/js/app.js` dosyasının başındaki `geometry()` ve `simulate()` fonksiyonlarındadır.
  Varsayılan roket `DEFAULTS`, hazır ayarlar `PRESETS` içindedir.
- Konu anlatımındaki bazı sayılar (örnek uçuşun puanı, ders kartlarındaki değerler) simülasyondan alınmıştır.
  Fizik modelini değiştirirseniz bu sayıları da güncelleyin.

## Not

Simülasyon karşılaştırma ve öğrenme amaçlı basitleştirilmiş bir modeldir; sonuçları gerçek testlerle doğrulayın.
Rehberdeki kural özetleri ASRW 2026 belgelerinin Türkçe çevirisine dayanır; resmî ve bağlayıcı olan
organizasyonun yayımladığı İngilizce belgelerdir.
