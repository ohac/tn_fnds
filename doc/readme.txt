UTAU用音声合成エンジン『tn_fnds』 version 0.0.9

■概要
このプログラムは森勢将雅氏のWORLD版UTAU合成エンジン『エターナルフォースブリサンプラー 
ジェントリー・ウィープス　〜相手は死ぬ，俺も死ぬ〜』(EFB-GW)をカスタマイズし、連続音や
子音速度に対応させたものです。
高速化やメモリ使用量削減、フラグ対応のために、合成エンジンのコア部分に手を加えていますが、
合成手法の基本的な考え方は変えていません。（変えていないつもりです。）

作成にあたっては、飴屋／菖蒲氏のworld4utauのソースコードを流用させていただいております。

■制約事項
・使用できないフラグが多数あります。
・合成処理の都合上、厳密なパラメータの解釈は行っていないため、合成された音声のタイミ
　ングやピッチが正確ではない場合があります。
・原音の状態により、合成後の音声に原音には無いノイズが入る場合があります。

■対応フラグ
・gフラグ
・tフラグ
・Bフラグ(BRE)
　50以上を設定したときに息成分が増えます。50未満を設定しても変化しません。
・bフラグ
　本家とは動作が違い、子音（無声部）を強調する動作を行います。
　b10で5%増、b100で500%増となります。
・Aフラグ(0~100 def0)
　独自フラグです。ピッチの変動に合わせて音量を変化させるます。ビブラート部分に
　ピンポイントで使用することをお勧めします。
・Oフラグ(-100~100 def0)
　独自フラグです。声の強さ（明るさ）を変更します。+方向を指定すると、低周波を抑制し
　高周波を増幅します。無声部へは影響しません。
・eフラグ
　独自フラグです。音符を引き伸ばす方法を変更します。パラメータはありません。
　デフォルトでは伸縮区間（原音設定の白い部分）をループさせて引き伸ばしますが、
　eフラグを指定すると伸縮区間をそのまま引き伸ばします（resamplerと同じ方式）。
　ループ方式の方が自然に聞こえることが多いですが、ループ特有のノイズが気になる
　場合はeフラグを使用してください。
・Wフラグ(-1, 50~1000)
　デスボイス用のフラグです。本家のWフラグとは動作が違います。
　F0分析の結果を指定した周波数で上書きします。-1を指定すると無声音として扱います。
　原音のプロフィールの『freq avg』の値を基本として調整してください。

■ライセンス
・GPLv3ライセンスとなっています。ソースコードを改変し配布する場合は、改変したソース
　コードも同様に配布してください。
　ライセンスの詳細はcopying.txtをご覧ください。

■更新履歴
2012/2/19  ver0.0.1 ・公開
                    ・EFB-GWからの変更点
                      合成に必要な範囲のみ分析するようにした
                      引き伸ばし範囲をUTAUの仕様にあわせた
                      ピッチベンドの適用方法をUTAUの仕様にあわせた
                      子音速度に対応した
                      その他、異常動作となる部分を複数修正した
2012/2/19  ver0.0.2 ・異常なピッチとなる不具合の修正
                    ・入力データによってフリーズすることがある不具合の修正
                    ・分析範囲を必要な範囲の前後100msに拡大した

2012/2/25  ver0.0.3 ・g,t,A,eフラグを実装
                    ・内部的なデータの持ち方を変えることにより高速化を図った
                    ・F0の平滑化、内部パラメータの調整により音質の向上を図った
                    ・ステレオの原音を読み込めるようにした

2012/3/10  ver0.0.4 ・Bフラグを実装
                    ・F0補正の改良による無声子音の高音質化、mod0での音程の改善
                    ・ピッチベンド適用方法をW4Uに合わせた（クロスフェード最適化が効きやすくなったかも？）

2012/3/17  ver0.0.5 ・b,Oフラグを実装
                    ・Bフラグのバランス調整

2012/3/31  ver0.0.6 ・Wフラグを実装

2012/6/9   ver0.0.7 ・内部処理の調整により音質の向上を図った
2012/6/9   ver0.0.8 ・内部処理の調整により音質の向上を図った（修正）
2013/4/11  ver0.0.9 ・FFTWを使わないようにした

■最新版入手先
　ZteerのUTAU関連物置き場
　http://z-server.game.coocan.jp/utau/utautop.html

■連絡先
　不具合などありましたら以下までご連絡ください。
・Zteer(BC_701@hotmail.com)
  Twitter(@zteer)

